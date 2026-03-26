package main

import (
	"crypto/sha256"
	"encoding/hex"
	"fmt"
	"os"
	"sort"
)

// hashSources computes a deterministic hash of source file paths and contents.
func hashSources(sources []string) (string, error) {
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
	return hex.EncodeToString(h.Sum(nil))[:16], nil
}
