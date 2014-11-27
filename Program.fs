open System
open System.Diagnostics
open Rosalind

let readia () = Console.ReadLine().Split() |> Array.map int

let lineSeq (r:System.IO.TextReader) :string seq =
        Seq.unfold(fun line -> 
            if line = null then 
                r.Close() 
                None 
            else 
                Some(line,r.ReadLine())) (r.ReadLine())

[<EntryPoint>]
let main argv =
(*
  // Test FASTA & GC
  let records = System.Console.In |> FASTA.parse
  let m = records |> Seq.maxBy (fun f -> f.dna |> DNA.gcCount)
  Console.WriteLine(sprintf "%s\n%f" (match m.id with | None -> "(none)" | Some i -> i) (m.dna |> DNA.gcCount))
*)

(*
  // Test protein translate
  let r = Console.ReadLine() |> RNA.translateToProtein
  match r with
    | (Some p, _) -> Console.WriteLine(p)
    | _           -> Console.WriteLine("Could not translate to a protein string")
*)

(*
  // Test DNA substrings
  let src=Console.ReadLine()
  let target=Console.ReadLine()
  let subs=DNA.find src target
  Console.WriteLine(sprintf "%A" subs)
*)

(*
  // Test Mortal Fibonacci Rabbits
  let [|n; m|]=readia()
  let a = DNA.mortalPopulation m |> Seq.nth (n-1)
  printfn "%u" a
*)

(*
  // Test Consensus
  let records = Console.In |> FASTA.parse |> Seq.cache
  let cm = records |> Seq.map (FASTA.dna) |> Consensus.calculate
  let profile = Consensus.profile cm

  printfn "%s" profile
  Consensus.print cm
*)

(*
  // RNA Splicing
  let rec removeIntrons (src:DNA.T) (introns:DNA.T list) =
    match introns with
      | intron::t as T -> let r = DNA.find src intron
                          match r with
                          | []    -> removeIntrons src t
                          | f::ft -> removeIntrons (src.Remove(f, String.length intron)) T
      | [] -> src
  
  let sin = new System.IO.StreamReader("/Users/marcusedwards/Projects/rosalind/data/rosalind_splc.txt")
  let records = sin |> FASTA.parse |> Seq.map (FASTA.dna) |> Seq.toList
  let extrons = removeIntrons (List.head records) (List.tail records)

  extrons
  |> RNA.transcribeDNA
  |> RNA.translateToProtein
  |> function | (Some p,_) -> p | (None, rna) -> sprintf "(Could not translate: %s)" rna
  |> printfn "%s"
*)

(*
  // Transition-Transversion ratio
  let sin = new System.IO.StreamReader(argv.[0])
      //if (Array.length argv > 0) then new System.IO.StreamReader(argv.[0])
      //else Console.In
  let records = sin |> FASTA.parse |> Seq.map (FASTA.dna) |> Seq.toArray

  let vd = DNA.visualDiff records.[0] records.[1]
  printfn "tr=%s\ntv=%s" (System.String.Concat (fst vd)) (System.String.Concat (snd vd))

  let d = DNA.diff records.[0] records.[1]
  printfn "tr=%u\ntv=%d" (fst d) (snd d)

  let trtv = DNA.trtvRatio records.[0] records.[1]
  printfn "%f" trtv
*)

(*
  // Test permutations
  let n = Console.ReadLine() |> int
  let ps = Util.permutations n |> Seq.cache

  printfn "%u" (Seq.length ps)
  ps |> Seq.iter (printfn "%A")
*)

(*
  // Test restriction sites: revp
  let sin = new System.IO.StreamReader(argv.[0])
  let dna = sin |> FASTA.parse |> Seq.map (FASTA.dna) |> Seq.head

  dna   
  |> DNA.findRestrictionSites
  |> Seq.iter (fun (i, site) -> printfn "%d %d" (i+1) (String.length site))
*)
 (* 
  // Test graph node degree: deg, ddeg
  printfn "argv=%A" argv
  let (src:System.IO.TextReader) = match Array.length argv with
                                   | 1 -> upcast new System.IO.StreamReader(argv.[0])
                                   | _ -> System.Console.In
  let header = src.ReadLine().Split();

  let graph =
    src
    |> lineSeq
    |> Seq.fold (fun g l -> let [|sn; e|] = l.Split() |> Array.map int
                            g |> Graph.AddUndirectedEdge sn e) Graph.Empty

  graph
  |> Graph.degrees
  |> List.iter (fun (_, count) -> printf "%u " count)
*)
(* Testing DNA.indexOf
   let src = DNA.fromString "GATTACA"

   Debug.Assert(Option.isNone (DNA.indexOf [||] src), "Failed: empty target")
   Debug.Assert(Option.isNone (DNA.indexOf src [||]), "Failed: empty source")

   let t3 = DNA.indexOf (DNA.fromString "GA") src
   Debug.Assert(Option.isSome t3 && Option.get t3 = 0, "Failed: didn't find leading characters");

   let t4 = DNA.indexOf (DNA.fromString "GATTACAGATTACA") src
   Debug.Assert(Option.isNone t4, "Failed: target is longer than source");

   let t5 = DNA.indexOf src src
   Debug.Assert(Option.isSome t5 && Option.get t5 = 0, "Failed: didn't find src in src");

   let t6 = DNA.indexOf (DNA.fromString "CA") src
   Debug.Assert(Option.isSome t6 && Option.get t6 = 5, "Failed: didn't find trailing characters");

   let t7 = DNA.indexOf (DNA.fromString "C") src
   Debug.Assert(Option.isSome t7 && Option.get t7 = 5, "Failed: didn't find single character");

   let t8 = DNA.indexOf (DNA.fromString "TTA") src
   Debug.Assert(Option.isSome t8 && Option.get t8 = 2, "Failed: didn't find  character in the middle");

   let t9 = DNA.indexOf (DNA.fromString "ACATTAG") src
   Debug.Assert(Option.isNone t9 , "Failed: full length non-matching string");

   let t10 = DNA.indexOf (DNA.fromString "XYZ") src
   Debug.Assert(Option.isNone t10, "Failed: non-matching string");
*)

   // Test longest substrings
   let (src:System.IO.TextReader) = match Array.length argv with
                                    | n when n >= 1 -> upcast new System.IO.StreamReader(argv.[0])
                                    | _             -> System.Console.In
//   let src = new System.IO.StreamReader("../../data/rosalind_lcsm.txt")
   let strings =
      src
      |> FASTA.parse
      |> Seq.map (FASTA.dna)
      |> List.ofSeq
      |> DNA.subStrings
   let lcsm =
      strings
      |> Seq.maxBy (DNA.length)
      |> DNA.toString


   printfn "longest sub-string=%s" lcsm

   0

