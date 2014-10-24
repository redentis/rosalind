module Rosalind.Consensus

type T = {A:int array;
          T:int array;
          G:int array;
          C:int array} with
    member this.length = Array.length this.A

let private createMatrix (length:int) :T =    
  {A=Array.create length 0;
   C=Array.create length 0;
   G=Array.create length 0;
   T=Array.create length 0}

let calculate (dnas:DNA.T seq) :T =
  let l = dnas |> Seq.head |> String.length
  let matrix = createMatrix l
  let incr (m:int array) i =
    m.[i] <- m.[i] + 1
  dnas
  |> Seq.iter (fun d -> d |> Seq.iteri (fun i bp -> match bp with
                                                    | 'A' -> incr matrix.A i
                                                    | 'C' -> incr matrix.C i
                                                    | 'G' -> incr matrix.G i
                                                    | 'T' -> incr matrix.T i
                                                    | _ -> ()))
  matrix
  
let print (cm:T) =
  let spacedString (cs:int array) = cs |> Array.fold (fun a i -> a + (sprintf "%u " i)) ""
  printfn "A: %s" (spacedString cm.A)
  printfn "C: %s" (spacedString cm.C)
  printfn "G: %s" (spacedString cm.G)
  printfn "T: %s" (spacedString cm.T)

let profile (cm:T) :DNA.T =
  [ for i in 0..(cm.length-1) do yield [('A',cm.A.[i]); ('C',cm.C.[i]); ('G',cm.G.[i]); ('T',cm.T.[i])]]
  |> List.fold (fun profile col -> let (bp, _) = col |> List.maxBy (fun (_,c)->c)
                                   bp::profile) []
  |> List.rev
  |> DNA.fromBasePairs

