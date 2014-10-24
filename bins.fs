open System

let readia () = Console.ReadLine().Split() |> Array.map int

[<EntryPoint>]
let main argv =
  let n=Console.ReadLine() |> int
  let m=Console.ReadLine() |> int
  let data=readia () |> Array.sort

  readia ()
  |> Array.iter (fun q -> match (Array.tryFindIndex (fun e->e=q) data) with
                          | Some i -> printf "%d " (i+1)
                          | None   -> printf "-1 ")
  0
