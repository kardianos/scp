package main

import (
	"crypto/sha256"
	"encoding/hex"
	"fmt"
	"os"
	"path/filepath"
	"regexp"
	"sort"
)

// hashBuild computes a deterministic hash of source file contents and the build command.
// Including the build command ensures that changing flags (e.g. -arch=sm_70 vs sm_89)
// invalidates the cache even when sources are unchanged.
func hashBuild(sources []string, cmd string) (string, error) {
	h := sha256.New()
	sorted := make([]string, len(sources))
	copy(sorted, sources)
	sort.Strings(sorted)

	for _, s := range sorted {
		data, err := os.ReadFile(s)
		if err != nil {
			return "", fmt.Errorf("read %s: %w", s, err)
		}
		h.Write([]byte(s))
		h.Write(data)
	}
	h.Write([]byte(cmd))
	return hex.EncodeToString(h.Sum(nil))[:16], nil
}

// localIncludeRe matches quoted local includes: #include "path/to/file.h".
var localIncludeRe = regexp.MustCompile(`(?m)^\s*#\s*include\s+"([^"]+)"`)

// sourceExts are file extensions scanned for local includes.
var sourceExts = map[string]bool{
	".c": true, ".cu": true, ".cc": true, ".cpp": true,
	".h": true, ".hh": true, ".hpp": true, ".cuh": true,
}

// expandSources scans the listed source files for local `#include "..."`
// directives (recursively) and appends any discovered files to the source set,
// so callers can pass just the main source file. Includes are resolved
// relative to the including file's directory first, then relative to the
// current working directory (repo layout). Originals keep their position;
// discovered files are appended, deduplicated.
func expandSources(sources []string) []string {
	seen := make(map[string]bool)
	var out []string
	add := func(p string) {
		p = filepath.Clean(p)
		if !seen[p] {
			seen[p] = true
			out = append(out, p)
		}
	}
	for _, s := range sources {
		add(s)
	}
	// Breadth-first over the growing list.
	for i := 0; i < len(out); i++ {
		src := out[i]
		if !sourceExts[filepath.Ext(src)] {
			continue
		}
		data, err := os.ReadFile(src)
		if err != nil {
			continue
		}
		dir := filepath.Dir(src)
		for _, m := range localIncludeRe.FindAllSubmatch(data, -1) {
			inc := string(m[1])
			if cand := filepath.Join(dir, inc); fileExists(cand) {
				add(cand)
			} else if fileExists(inc) {
				add(inc)
			}
		}
	}
	return out
}

// compileSources filters a source list down to compilable files (.c/.cu/.cc/.cpp),
// excluding headers that are present only for upload/hashing.
func compileSources(sources []string) []string {
	var out []string
	for _, s := range sources {
		switch filepath.Ext(s) {
		case ".c", ".cu", ".cc", ".cpp":
			out = append(out, s)
		}
	}
	return out
}

// fileExists reports whether path exists and is a regular file.
func fileExists(path string) bool {
	st, err := os.Stat(path)
	return err == nil && st.Mode().IsRegular()
}
