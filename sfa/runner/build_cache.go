package main

import (
	"crypto/sha256"
	"encoding/hex"
	"fmt"
	"os"
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
