module Rosalind.FASTA

  open System.IO

  type T = {id:string option; description:string option; dna:DNA.T}

  let id (e:T) = e.id

  let description (e:T) = e.description

  let dna (e:T) = e.dna

  let lineSeq (r:TextReader) :string seq =
        Seq.unfold(fun line -> 
            if line = null then 
                r.Close() 
                None 
            else 
                Some(line,r.ReadLine())) (r.ReadLine())
                
  let parse (r:TextReader) :T seq =
    let parseHeader () =
      let line=r.ReadLine();
      match line.IndexOf(' ') with
        | -1 -> Some line.[1..],  None
        | n  -> Some line.[1..n], Some line.[n+1..]
    let parseSequences () =
      let rec aux (dnas:string list) =
        match r.Peek() with
          | -1 | 62 -> dnas |> List.rev |> List.reduce (+)
          | _ -> let dna=r.ReadLine()
                 aux (dna::dnas)
      aux []  
    let rec scanLines () =
        seq {           
            match r.Peek() with
            | -1  -> ()
            | 62  -> let i, d = parseHeader ()
                     let dna = parseSequences ()
                     yield {id=i; description=d; dna=DNA.fromString dna}
                     yield! scanLines ()
            | _   -> r.ReadLine() |> ignore
                     yield! scanLines ()
        }
    scanLines ()
