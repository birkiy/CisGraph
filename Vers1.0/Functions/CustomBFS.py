
from Functions.Packages import *


def generic_bfs_edgesCustom(G, source, neighbors=None, depth_limit=None, currentDepth=None, distance=None, shell=None):
    visited = {source}
    if depth_limit is None:
        depth_limit = len(G)
    if distance == False:
        queue = deque([(0, source, depth_limit, neighbors(source))])
    else:
        queue = [(0, source, depth_limit, neighbors(source))]
        heapq.heapify(queue)
    while queue:
        distanceR, parent, depth_now, children = queue[0]
        try:
            child = next(children)
            currentDepth += [depth_now]
            if child not in visited:
            # if child not in visited and G.nodes[child]["nodeClass"] != "upP" and G.nodes[child]["nodeClass"] != "dwP" and G.nodes[child]["nodeClass"] != "hom":
                shell[0] = [source]
                if shell[depth_limit - depth_now + 1] is None:
                    shell[depth_limit - depth_now + 1] = [child]
                else:
                    shell[depth_limit - depth_now + 1] += [child]
                yield parent, child
                visited.add(child)
                if depth_now > 1:
                    if distance == False:
                        queue.append((distanceR, child, depth_now - 1, neighbors(child)))
                    else:
                        if G.edges()[(parent, child)]["distance"] == "INF":
                            distanceNow = 9999999999999999999999
                        else:
                            distanceNow = abs(float(G.edges()[(parent, child)]["distance"]))
                        heapq.heappush(queue, (distanceR + distanceNow, child, depth_now - 1, neighbors(child)))
        except StopIteration:
            if distance == False:
                queue.popleft()
            else:
                heapq.heappop(queue)




def bfs_edgesCustom(G, source, reverse=False, depth_limit=None, current_depth=None, distance=None, shell=None):
    if reverse and G.is_directed():
        successors = G.predecessors
    else:
        successors = G.neighbors
    for e in generic_bfs_edgesCustom(G, source, successors, depth_limit, currentDepth=current_depth, distance=distance, shell=shell):
        yield e


def bfs_treeCustom(G, source, reverse=False, depth_limit=None, current_depth=None, distance=None, shell=None):
    T = nx.DiGraph()
    T.add_node(source)
    edges_gen = bfs_edgesCustom(G, source, reverse=reverse, depth_limit=depth_limit, current_depth=current_depth, distance=distance, shell=shell)
    T.add_edges_from(edges_gen)
    return T
