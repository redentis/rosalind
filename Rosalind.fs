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

module DNA = 
   open System
   
   // Base type for representing DNA and use in DNA operations
   type T = string

   type codon =
       | T
       | A
       | C
       | G
       | Unknown of char
   
   type T' = codon seq

   type Mutation = 
      | Transition of char
      | Transversion of (char * char)
   
   let fromString (s : string) : T' = 
      s.ToCharArray()
      |> Seq.map (function | 'T' -> T | 'A' -> A | 'G' -> G | 'C' -> C | _ as c -> Unknown c) 

   let fromBasePairs (bps : char list) : T = System.String.Concat bps

   let length (dna: T) :int =
       String.length dna
   
   // Return the reverse compliment of the given DNA string.
   let reverseComplement (dna : T) : T = 
      dna
      |> Seq.toList
      |> List.rev
      |> List.map (function 
            | 'T' -> 'A'
            | 'A' -> 'T'
            | 'C' -> 'G'
            | 'G' -> 'C'
            | b -> b)
      |> System.String.Concat
   
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
         |> Seq.filter (fun bp -> bp = 'G' || bp = 'C')
         |> length
      
      let l = dna |> Seq.length
      (float gc) / (float l) * 100.00
   
   // Search for the given target sequence in the source sequence and return the indexes in
   // the DNA string where it occurs.
   let find (src : T) (target : T) : int list = 
      let l = (String.length src) - 1
      let n = (String.length target) - 1
      
      let rec aux i is = 
         if i + n > l then is |> List.rev
         else if src.[i..(i + n)] = target then aux (i + 1) (i :: is)
         else aux (i + 1) is
      aux 0 []
   
   let remove (src : T) (target : T) : T = 
      let i = find src target |> List.head
      src.Remove(i, target.Length)
   
   // Calculate the difference between two DNA sequences returning a tuple with
   // transition count and transversion count. The two values combined give the
   // hamming distance between the two sequences.
   let diff (l : T) (r : T) : int * int = 
      Seq.zip l r |> Seq.fold (fun (tr, tv) p -> 
                        match p with
                        | ('A', 'G') | ('G', 'A') | ('C', 'T') | ('T', 'C') -> (tr + 1, tv)
                        | (a, b) when a = b -> (tr, tv)
                        | _ -> (tr, tv + 1)) (0, 0)
   
   let visualDiff (l : T) (r : T) : char list * char list = 
      let (tr, tv) = 
         Seq.zip l r |> Seq.fold (fun (tr, tv) p -> 
                           match p with
                           | ('A', 'G') | ('G', 'A') | ('C', 'T') | ('T', 'C') -> ('.' :: tr, ' ' :: tv)
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
                  let candidate = left + right
                  // printfn "candidate=%A" candidate
                  if (dna.[i..bound] = candidate) then yield (i, candidate)
      }

module Protein = 
   type T = string

module RNA = 
   type T = string
   
   type Codon = 
      | AminoAcid of char
      | Stop
   
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
      dna |> String.map (function 
                | 'T' -> 'U'
                | b -> b)
   
   let translateToProtein (rna : T) : Protein.T option * T = 
      let l = String.length rna
      
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
            match codonMap.TryFind codon with
            | None | Some Stop -> attach (j + 1)
            | Some(AminoAcid a) -> translate (a :: alist) (j + 1)
      translate [] 0
