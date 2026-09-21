// Package qgraph implements the quantized harmonic graph (QHG): a discrete
// law in which emission, absorption, bound shells and charge are the same
// mechanism read at different moments of one ladder.
//
// Ontology (v93 constraints):
//   - one law: energy changes mode, never vanishes;
//   - no background: the graph is the CONTACT graph of a live roster, it is
//     re-derivable at every instant and carries no immortal index set;
//   - no imported U(1) / electron species: charge is a spectral occupancy
//     defect of the graph's own harmonics;
//   - free-cell is the physical substrate; this package is the symbolic one.
//
// See docs/UNIFIED_LAW.md.
package qgraph

// Graph is a finite undirected weighted simple graph.
//
// It is a chart of a contact set, not a lattice: nothing here assumes the
// vertex set survives to the next instant. Every quantity the law uses is a
// spectral or topological invariant of the graph as it stands NOW.
type Graph struct {
	N     int         // vertices
	Adj   [][]int     // adjacency list
	W     [][]float64 // symmetric weights; 0 means "no edge"
	Edges [][2]int    // undirected edge list, i<j
}

// NewGraph returns an empty graph on n vertices.
func NewGraph(n int) *Graph {
	g := &Graph{N: n, Adj: make([][]int, n), W: make([][]float64, n)}
	for i := range g.W {
		g.W[i] = make([]float64, n)
	}
	return g
}

// AddEdge adds an undirected edge of weight w. Repeat edges are ignored.
func (g *Graph) AddEdge(i, j int, w float64) {
	if i == j || i < 0 || j < 0 || i >= g.N || j >= g.N {
		return
	}
	if i > j {
		i, j = j, i
	}
	if g.W[i][j] != 0 {
		return
	}
	g.Edges = append(g.Edges, [2]int{i, j})
	g.Adj[i] = append(g.Adj[i], j)
	g.Adj[j] = append(g.Adj[j], i)
	g.W[i][j] = w
	g.W[j][i] = w
}

// NewPath builds the path P_n (an open chain).
func NewPath(n int) *Graph {
	g := NewGraph(n)
	for i := 0; i < n-1; i++ {
		g.AddEdge(i, i+1, 1)
	}
	return g
}

// NewCycle builds the ring C_n. This is the programme's own retention object
// (v93 RING_DNLS: the C6 ring that held W=+2 in vacuum).
func NewCycle(n int) *Graph {
	g := NewPath(n)
	if n >= 3 {
		g.AddEdge(0, n-1, 1)
	}
	return g
}

// NewComplete builds K_n.
func NewComplete(n int) *Graph {
	g := NewGraph(n)
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			g.AddEdge(i, j, 1)
		}
	}
	return g
}

// NewHypercube builds Q_d, the d-dimensional hypercube graph (2^d vertices).
// Its Laplacian spectrum is 2k with multiplicity C(d,k): a clean shell ladder.
func NewHypercube(d int) *Graph {
	n := 1 << d
	g := NewGraph(n)
	for v := 0; v < n; v++ {
		for b := 0; b < d; b++ {
			g.AddEdge(v, v^(1<<b), 1)
		}
	}
	return g
}

// NewTorus builds the discrete d-torus with side L (L^d vertices, periodic).
// Its spectral dimension is d, which is what fixes the Stefan exponent.
func NewTorus(d, L int) *Graph {
	n := 1
	for i := 0; i < d; i++ {
		n *= L
	}
	g := NewGraph(n)
	stride := make([]int, d)
	s := 1
	for i := 0; i < d; i++ {
		stride[i] = s
		s *= L
	}
	for v := 0; v < n; v++ {
		for i := 0; i < d; i++ {
			c := (v / stride[i]) % L
			u := v + ((c+1)%L-c)*stride[i]
			g.AddEdge(v, u, 1)
		}
	}
	return g
}

// NewIcosahedron builds the icosahedral graph (12 vertices, 5-regular).
// Its Laplacian degeneracies are 1,3,5,3 — an s/p/d-like ladder produced by
// graph symmetry alone, with no orbital imported.
func NewIcosahedron() *Graph {
	// Standard 12-vertex, 30-edge adjacency (vertex 0 and 11 are the poles).
	adj := [][]int{
		{1, 2, 3, 4, 5},
		{0, 2, 5, 6, 7},
		{0, 1, 3, 7, 8},
		{0, 2, 4, 8, 9},
		{0, 3, 5, 9, 10},
		{0, 1, 4, 6, 10},
		{1, 5, 7, 10, 11},
		{1, 2, 6, 8, 11},
		{2, 3, 7, 9, 11},
		{3, 4, 8, 10, 11},
		{4, 5, 6, 9, 11},
		{6, 7, 8, 9, 10},
	}
	g := NewGraph(12)
	for i, row := range adj {
		for _, j := range row {
			g.AddEdge(i, j, 1)
		}
	}
	return g
}

// Laplacian returns the combinatorial Laplacian L = D - W.
//
// L is the discrete d†d of the graph's own exterior derivative: it is the
// only second-order operator the contact set defines without extra structure.
// Its eigenvectors are THE harmonics of this note.
func (g *Graph) Laplacian() [][]float64 {
	L := make([][]float64, g.N)
	for i := 0; i < g.N; i++ {
		L[i] = make([]float64, g.N)
		deg := 0.0
		for j := 0; j < g.N; j++ {
			if g.W[i][j] != 0 {
				L[i][j] = -g.W[i][j]
				deg += g.W[i][j]
			}
		}
		L[i][i] = deg
	}
	return L
}

// CycleRank is the first Betti number b1 = |E| - |V| + components.
// It is the dimension of the space in which winding (and therefore charge)
// can live: a forest has b1 = 0 and can carry no circulation at all.
func (g *Graph) CycleRank() int {
	return len(g.Edges) - g.N + g.Components()
}

// Components counts connected components.
func (g *Graph) Components() int {
	seen := make([]bool, g.N)
	c := 0
	for s := 0; s < g.N; s++ {
		if seen[s] {
			continue
		}
		c++
		stack := []int{s}
		seen[s] = true
		for len(stack) > 0 {
			v := stack[len(stack)-1]
			stack = stack[:len(stack)-1]
			for _, u := range g.Adj[v] {
				if !seen[u] {
					seen[u] = true
					stack = append(stack, u)
				}
			}
		}
	}
	return c
}

// IsForest reports whether the graph has no cycles (b1 == 0).
//
// Theorem (radiation neutrality, UNIFIED_LAW §6.3): if the exterior sector a
// quantum escapes into is a forest, it can carry no winding, so no emission
// or absorption event can change the circulation books of the emitter.
func (g *Graph) IsForest() bool { return g.CycleRank() == 0 }
