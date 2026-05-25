import ShulgaParameters

def main : IO Unit := do
  let lambda := SCPv59.Furey7D.Shulga.lambda_raw_exact
  let mu := SCPv59.Furey7D.Shulga.mu_raw_exact
  let ratio := SCPv59.Furey7D.Shulga.geometric_ratio
  
  -- Create a basic JSON string output
  let jsonStr := 
    s!"\{\n" ++
    s!"  \"lambda_raw_exact\": \"{lambda.num}/{lambda.den}\",\n" ++
    s!"  \"mu_raw_exact\": \"{mu.num}/{mu.den}\",\n" ++
    s!"  \"geometric_ratio\": \"{ratio.num}/{ratio.den}\"\n" ++
    s!"}"
  
  IO.println jsonStr
