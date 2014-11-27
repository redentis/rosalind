namespace Rosalind

module Util = 
   let reverseJoinF xs f = 
      xs
      |> List.rev
      |> List.reduce f
   
   let reverseJoin xs = reverseJoinF xs (+)
   
   let intersperse c cs = 
      let rec aux = 
         function 
         | [] -> []
         | [ h; t ] -> h :: c :: t :: []
         | h :: t -> h :: c :: (aux t)
      aux cs
   
   let rec private insertions x = 
      function 
      | [] -> [ [ x ] ]
      | (y :: ys) as l -> (x :: l) :: (List.map (fun x -> y :: x) (insertions x ys))
   
   let rec private perms = 
      function 
      | [] -> seq [ [] ]
      | x :: xs -> Seq.concat (Seq.map (insertions x) (perms xs))
   
   let permutations (n : int) : int list seq = perms [ 1..n ]

   let rec public take n = function
      | [] -> []
      | h::t as _xs -> if n=0 then [] else h::(take (n-1) t)

   let rec public drop n = function
      | [] -> []
      | h::t as _xs -> if n=0 then _xs else drop (n-1) t

   let public remove i n xs = (take i xs)@(drop (i+n)xs)

module DNA = 
   open System

   type DNABase =
       | T
       | A
       | C
       | G
       | Unknown of char

   // Base type for representing DNA and use in DNA operations
   type T = DNABase array

   type Mutation = 
      | Transition of char
      | Transversion of (char * char)
     

   let toString (dna:T) :string = 
       dna |> Array.map (function | T -> 'T' | A -> 'A' | G -> 'G' | C -> 'C' | Unknown c -> c)
           |> System.String.Concat

   let fromBasePairs (bps:char array) :T =
      let dna = Array.create bps.Length (DNABase.Unknown '.')
      bps
      |> Seq.iteri (fun i -> function
                              | 'T' -> dna.[i] <- T
                              | 'A' -> dna.[i] <- A
                              | 'G' -> dna.[i] <- G
                              | 'C' -> dna.[i] <- C
                              | _ as c -> dna.[i] <- Unknown c)
      dna

   // Creates a new DNA.T from the given string
   let fromString (s : string) : T = 
      s.ToCharArray() |> fromBasePairs

   // Returns the length of the DNA string 
   let length (dna: T) :int =
       dna.Length
   
   // Searches src for the first occurance of target. If an occurance is found, returns
   // an int option with the index of the first character.
   let indexOf (target:T) (src:T) :int option =
      match target, src with 
      | [||], _ | _, [||] -> None
      | _                 -> let mt = (length target)-1
                             let ms = (length src)-1
                             let rec aux t s :int option =
                                if t > mt then Some (s-mt-1)
                                else if s > ms then None
                                else
                                  if target.[t] = src.[s] then aux (t+1) (s+1) else aux 0 (s+1)
                             aux 0 0

   // Return the reverse compliment of the given DNA string.
   let reverseComplement (dna : T) : T = 
      dna
      |> Array.rev
      |> Array.map (function 
                     | T -> A
                     | A -> T
                     | C -> G
                     | G -> C
                     | Unknown _ as c -> c)
   
   // Calculate the probability of a random offspring from a population of
   // k homozygous dominant, m heterozygous, n homozygous recessive expressing
   // a dominant allele.
   let pOffspringDominant (k : int) (m : int) (n : int) : float = 
      let Total = k + m + n
                  |> float
      let Tp = Total ** 2.0 - Total
      let fk = float k
      let fm = float m
      let fn = float n
      (fk ** 2.0 - fk) / Tp + (2.0 * fk * fm) / Tp + (2.0 * fk * fn) / Tp + (fm * fn) / Tp 
      + ((fm ** 2.0 - fm) / Tp) * 0.75
   
   // Generate a sequence representing the population size in each period of reproduction
   // assuming that each reproductive pair has k offspring in each generation.
   let population (growth : int) : uint64 seq = 
      let growthp = uint64 growth
      Seq.unfold (fun (infant, mature) -> 
         let births = mature * growthp
         Some((infant + mature), (births, infant + mature))) (1UL, 0UL)
      |> Seq.cache
   
   // Generate a sequence repsenting the populate size in each period of reproduction assuming
   // that each reproductive pair has a given lifespan. This method assumes that each pair
   // only produces a single reproductive pair in each period.
   let mortalPopulation (lifespan : int) : uint64 seq = 
      let generations = Array.create lifespan 0UL
      Seq.unfold (fun (infant : uint64, mature : uint64, gs : uint64 array, i) -> 
         //                  Console.WriteLine(sprintf "infant=%u; mature=%u; gs=%A; i=%u" infant mature gs i)
         let dead = gs.[i]
         let survivors = mature - dead
         let ii = (i + 1) % lifespan
         gs.[i] <- infant
         Some((survivors + infant), (survivors, (survivors + infant), gs, ii))) (1UL, 0UL, generations, 0)
      |> Seq.cache
   
   // Calculate the GC count for a given DNA string as a percentage between 0 and 100.
   let gcCount (dna : T) = 
      let gc = 
         dna
         |> Array.filter (fun bp -> bp = G || bp = C)
         |> Array.length
      
      (float gc) / (float dna.Length) * 100.00
   
   // Search for the given target sequence in the source sequence and return the indexes in
   // the DNA string where it occurs.
   let find (src : T) (target : T) : int list = 
  
      let l = src.Length - 1
      let n = target.Length - 1
      
      let rec aux i is = 
         if i + n > l then is |> List.rev
         else if src.[i..(i + n)] = target then aux (i + 1) (i :: is)
         else aux (i + 1) is
      aux 0 []
   
   let remove (src : T) (target : T) : T = 
      match (find src target) with
      | []   -> src
      | h::_ -> Array.append src.[..h] src.[(h+target.Length)..]
   
   // Calculate the difference between two DNA sequences returning a tuple with
   // transition count and transversion count. The two values combined give the
   // hamming distance between the two sequences.
   let diff (l : T) (r : T) : int * int = 
      Seq.zip l r |> Seq.fold (fun (tr, tv) p -> 
                        match p with
                        | (A, G) | (G, A) | (C, T) | (T, C) -> (tr + 1, tv)
                        | (a, b) when a = b -> (tr, tv)
                        | _ -> (tr, tv + 1)) (0, 0)
   
   let visualDiff (l : T) (r : T) : char list * char list = 
      let (tr, tv) = 
         Seq.zip l r |> Seq.fold (fun (tr, tv) p -> 
                           match p with
                           | (A, G) | (G, A) | (C, T) | (T, C) -> ('.' :: tr, ' ' :: tv)
                           | (a, b) when a = b -> (' ' :: tr, ' ' :: tv)
                           | _ -> (' ' :: tr, '.' :: tv)) ([], [])
      ((List.rev tr), (List.rev tv))
   
   // Calculate the error between two DNA strings l and r
   let hamm (l : T) (r : T) : int = 
      let (tr, tv) = diff l r
      tr + tv
   
   // Calculate the transition-transversion ratio of two DNA sequences
   let trtvRatio (l : T) (r : T) : float = 
      let (tr, tv) = diff l r
      (float tr) / (float tv)
   
   let findRestrictionSites (dna : T) : (int * T) seq = 
      let max = (Seq.length dna) - 1
      seq { 
         for i in 0..max do
            for l in 2..6 do
               let mirror = i + l - 1
               let bound = mirror + l
               // printfn "i=%d; l=%d; mirror=%d; bound=%d" i l mirror bound
               if (mirror <= max) && (bound <= max) then 
                  let left = dna.[i..mirror]
                  let right = reverseComplement left 
                  let candidate = Array.append left right 
                  // printfn "candidate=%A" candidate
                  if (dna.[i..bound] = candidate) then yield (i, candidate)
      }

   let subStrings (dnas:T list) =
      let initialCandidates (a:T) (b:T) :T list =
         printfn "a=%s; b=%s" (toString a) (toString b)
         let l = length b - 1
         [ for i in 0..(l-1) do
              for j in (i+1)..l do
                 let test = b.[i..j]
                 printfn "   testing %s" (toString test)
                 if Option.isSome (indexOf test a) then yield test]
      let rec search (_dnas:T list) (candidates:T list) :T list =
         match _dnas, candidates with
         | _, []   -> []
         | [], c   -> c
         | h::t, c -> search t (c |> List.filter (fun e -> match (indexOf e h) with | Some _ -> true | None -> false))
      match dnas with
      | [] | [_]-> []
      | a::b::r -> let c = (initialCandidates a b)@(initialCandidates b a)
                   c |> List.iter (fun e -> printfn "%s" (toString e))
                   search r c
module Protein = 
   type T = string

module RNA = 

   type RNABase =
      | G
      | U
      | A
      | C
      | Unknown of char

   type T = RNABase array          
   
   type Codon = 
      | AminoAcid of char
      | Stop

   let fromString (s : string) : T = 
      let rna = Array.create s.Length (RNABase.Unknown '.')
      s
      |> Seq.iteri (fun i -> function
                              | 'U' -> rna.[i] <- U
                              | 'A' -> rna.[i] <- A
                              | 'G' -> rna.[i] <- G
                              | 'C' -> rna.[i] <- C
                              | _ as c -> rna.[i] <- Unknown c)
      rna

   let toString (rna:T) :string = 
       rna |> Array.map (function | U -> 'U' | A -> 'A' | G -> 'G' | C -> 'C' | Unknown c -> c)
           |> System.String.Concat

   let length (rna:T) :int =
      rna.Length

   let codonMap = 
      Map.empty.Add("UUU", AminoAcid 'F').Add("UUC", AminoAcid 'F').Add("UUA", AminoAcid 'L').Add("UUG", AminoAcid 'L')
         .Add("UCU", AminoAcid 'S').Add("UCC", AminoAcid 'S').Add("UCA", AminoAcid 'S').Add("UCG", AminoAcid 'S')
         .Add("UAU", AminoAcid 'Y').Add("UAC", AminoAcid 'Y').Add("UAA", Stop).Add("UAG", Stop)
         .Add("UGU", AminoAcid 'C').Add("UGC", AminoAcid 'C').Add("UGA", Stop).Add("UGG", AminoAcid 'W')
         .Add("CUU", AminoAcid 'L').Add("CUC", AminoAcid 'L').Add("CUA", AminoAcid 'L').Add("CUG", AminoAcid 'L')
         .Add("CCU", AminoAcid 'P').Add("CCC", AminoAcid 'P').Add("CCA", AminoAcid 'P').Add("CCG", AminoAcid 'P')
         .Add("CAU", AminoAcid 'H').Add("CAC", AminoAcid 'H').Add("CAA", AminoAcid 'Q').Add("CAG", AminoAcid 'Q')
         .Add("CGU", AminoAcid 'R').Add("CGC", AminoAcid 'R').Add("CGA", AminoAcid 'R').Add("CGG", AminoAcid 'R')
         .Add("AUU", AminoAcid 'I').Add("AUC", AminoAcid 'I').Add("AUA", AminoAcid 'I').Add("AUG", AminoAcid 'M')
         .Add("ACU", AminoAcid 'T').Add("ACC", AminoAcid 'T').Add("ACA", AminoAcid 'T').Add("ACG", AminoAcid 'T')
         .Add("AAU", AminoAcid 'N').Add("AAC", AminoAcid 'N').Add("AAA", AminoAcid 'K').Add("AAG", AminoAcid 'K')
         .Add("AGU", AminoAcid 'S').Add("AGC", AminoAcid 'S').Add("AGA", AminoAcid 'R').Add("AGG", AminoAcid 'R')
         .Add("GUU", AminoAcid 'V').Add("GUC", AminoAcid 'V').Add("GUA", AminoAcid 'V').Add("GUG", AminoAcid 'V')
         .Add("GCU", AminoAcid 'A').Add("GCC", AminoAcid 'A').Add("GCA", AminoAcid 'A').Add("GCG", AminoAcid 'A')
         .Add("GAU", AminoAcid 'D').Add("GAC", AminoAcid 'D').Add("GAA", AminoAcid 'E').Add("GAG", AminoAcid 'E')
         .Add("GGU", AminoAcid 'G').Add("GGC", AminoAcid 'G').Add("GGA", AminoAcid 'G').Add("GGG", AminoAcid 'G')
   
   let transcribeDNA (dna : DNA.T) : T = 
      dna |> Array.map (function 
                | DNA.DNABase.T -> U
                | DNA.DNABase.G -> G
                | DNA.DNABase.A -> A
                | DNA.DNABase.C -> C
                | DNA.DNABase.Unknown c -> Unknown c)
   
   let translateToProtein (rna : T) : Protein.T option * T = 
      let l = rna.Length
      
      let rec translate (alist : char list) i = 
         let attach t = 
            match alist with
            | [] -> (None, rna.[t..])
            | _ -> 
               (Some(alist
                     |> List.rev
                     |> System.String.Concat), rna.[t..])
         match i with
         | i when i >= l -> attach l
         | i -> 
            let j = i + 2
            let codon = rna.[i..j]
            match codonMap.TryFind (codon |> toString) with
            | None | Some Stop -> attach (j + 1)
            | Some(AminoAcid a) -> translate (a :: alist) (j + 1)
      translate [] 0
