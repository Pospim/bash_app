#!/usr/bin/env python3


class GraphCleaningDFSException(Exception):
    def __init__(self, DFSnode, BFSnode):
        super().__init__(f"DFS from {DFSnode} encountered a node not visited during BFS from {BFSnode}.")
        self.DFSnode = DFSnode
        self.BFSnode = BFSnode
