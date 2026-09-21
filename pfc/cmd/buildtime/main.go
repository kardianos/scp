
package main
import (
  "fmt"
  "time"
  "scp/pfc/core"
)
func main() {
  law := core.DefaultLaw()
  law.L = 120
  t0 := time.Now()
  w := core.BuildLattice(law, 1728000)
  fmt.Printf("NC=%d slots=%d build=%v\n", w.NC(), len(w.Slots), time.Since(t0))
  t1 := time.Now()
  for i:=0;i<3;i++ { w.Step() }
  fmt.Printf("3 steps=%v\n", time.Since(t1))
}
