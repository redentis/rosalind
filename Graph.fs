module Rosalind.Graph

// Edge (target, weight)
type 'U edge = 
    | Edge of 'U
    | WeightedEdge of 'U * int

// Node (id, edge-list)
type 'U node = Node of 'U * edge<'U> list

type 'U T = T of node<'U> list

let nodeId = function 
             | Node (id, _) -> id

let edges  = function
             | Node (_, edges)-> edges

// ---------------------------------

let Empty:T<_> = (T[])

let AddNode<'U when 'U : equality> node (T(graph)) : T<'U> = 
    let (Node(nid, _)) = node
    graph
    |> List.tryFind (fun (Node(id, _)) -> id=nid)
    |> function
       | None   -> T(node::graph)
       | Some _ -> T(graph)

let RemoveNode (Node(id,_)) graph :T<'a> = 
    let rec aux ns = 
        match ns with
        | [] -> []
        | h::t when nodeId h = id -> t
        | h::t -> h::aux t
    let (T(nodes)) = graph
    nodes |> aux |> T

let AddDirectedEdge (srcId:'a) (targetId:'a) (graph:T<'a>) : T<'a> =
    let (T(nodes)) = graph
    let _graph = 
        nodes
        |> List.tryFind (fun (Node(id,_)) -> id = srcId)
        |> function 
            | None                        -> graph |> AddNode (Node(srcId, [Edge targetId ]))
            | Some (Node(id, edges) as n) -> graph |> RemoveNode n |> AddNode (Node (id, (Edge targetId)::edges))    
    nodes
    |> List.tryFind (fun (Node(id, _)) -> id = targetId)
    |> function
       | None  -> _graph |> AddNode (Node(targetId, []))
       | _     -> _graph
   
let AddUndirectedEdge<'a when 'a : equality> (srcId:'a) (targetId:'a) (graph:T<'a>) : T<'a> =
    graph
    |> AddDirectedEdge srcId targetId
    |> AddDirectedEdge targetId srcId

let degrees (T(nodes)) :('a * int) list =
    nodes
    |> List.map (fun (Node(id, edges)) -> (id, List.length edges))
    |> List.sort
   