



from Functions.Packages import *


def dfs_edgesCustom(G, source=None, depth_limit=None, current_depth=None):

    if source is None:
        # edges for all components
        nodes = G
    else:
        # edges for components with source
        nodes = [source]
    visited = set()
    if depth_limit is None:
        depth_limit = len(G)
    for start in nodes:
        if start in visited:
            continue
        visited.add(start)
        stack = [(start, depth_limit, iter(G[start]))]
        while stack:
            parent, depth_now, children = stack[-1]
            try:
                child = next(children)
                # current_depth += [depth_now]
                if child not in visited and G.nodes[child]["nodeClass"] != "dwP" and G.nodes[child]["nodeClass"] != "hom":
                    yield parent, child
                    visited.add(child)
                    if depth_now > 1:
                        stack.append((child, depth_now - 1, iter(G[child])))
                elif G.nodes[child]["nodeClass"] == "upP":
                    return visited

            except StopIteration:
                stack.pop()
        break
