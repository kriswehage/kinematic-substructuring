#!/usr/bin/env python3

'''
kinematic_substructuring.py

This script demonstrates the kinematic substructuring algorithm
described in Generalized Coordinate Partitioning for Complex
Mechanisms Based on Kinematic Substructuring.

Copyright(c) 2015 Kristopher Wehage
All rights reserved.

Redistribution and use in source and binary forms, with or
without modification, are permitted provided that the following
conditions are met:

1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above
copyright notice, this list of conditions and the following
disclaimer in the documentation and/or other materials provided
with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
SUCH DAMAGE.
'''

__author__ = '''Kristopher Wehage (ktwehage@ucdavis.edu)'''


def get_key(item):
    return item[0]


class SpanningTree(object):
    def __init__(self, P, loop_terminators, parent_joints):
        self.parent_joints = parent_joints
        self.P = P
        self.I = range(len(P))
        self.loop_terminators = loop_terminators
        self.L = []
        self.R = []
        self.loop_base_nodes = []
        self.levels = []
        self.substructures = []

        self.generate_left_right_list()
        self.generate_level_list()
        self.generate_loop_member_list()
        self.merge_loops()
        self.reorder_bodies_joints()

    def generate_left_right_list(self):
        # generate L & R from P.
        self.L = [0] * len(self.P)
        self.R = [0] * len(self.P)
        for j in range(0, len(self.I)):
            i = self.P[j]
            if self.L[i] == 0:
                self.L[i] = j
            else:
                k = self.L[i]
                while self.R[k] > 0:
                    k = self.R[k]
                self.R[k] = j

    def in_path(self, j, k):
        # is k in the kinematic path of j?
        in_path = False
        while j > 0:
            if j is k:
                in_path = True
                break
            j = P[j]
        return in_path

    def generate_level_list(self):
        # generate level list using preorder tree traversal
        S = []
        self.levels = [0] * len(self.I)
        k = self.L[0]
        while(S or k > 0):
            if(k > 0):
                self.levels[k] = self.levels[self.P[k]] + 1
                if (self.R[k] > 0):
                    S.append(self.R[k])
                k = self.L[k]
            else:
                k = S.pop()

    def generate_loop_member_list(self):
        # list generation and sorting are both performed in this method
        self.loop_members = [[] for i in range(len(self.loop_terminators))]
        for i in range(len(self.loop_terminators)):
            positive = self.loop_terminators[i][0]
            negative = self.loop_terminators[i][1]

            self.loop_members[i].append(negative)
            self.loop_members[i].append(positive)

            while self.levels[positive] > self.levels[negative]:
                positive = self.P[positive]
                self.loop_members[i].append(positive)

            while self.levels[negative] > self.levels[positive]:
                negative = self.P[negative]
                self.loop_members[i].append(negative)

            while positive is not negative:
                negative = self.P[negative]
                positive = self.P[positive]
                self.loop_members[i].append(negative)
                self.loop_members[i].append(positive)

            self.loop_members[i] = self.loop_members[i][0:-1][::-1]
        self.loop_members = sorted(self.loop_members, key=get_key)

    def merge_loops(self):
        self.substructures.append([])
        self.substructures[0].append(self.loop_members.pop(0))
        while self.loop_members:
            mergeCandidate = self.loop_members.pop(0)
            merge = False
            for substructure in self.substructures:
                for loop in substructure:
                    if self.in_path(mergeCandidate[0], loop[0]):
                        if any([e in loop[1:] for e in mergeCandidate[1:]]):
                            substructure.append(mergeCandidate)
                            merge = True
                            break
            if not merge:
                S = []
                S.append(mergeCandidate)
                self.substructures.append(S)

    def reorder_bodies_joints(self):
        # flatten substructures and sort members of each kinematic substructure
        # by level in spanning tree
        flattened_substructures = []
        for substructure in self.substructures:
            flattened_substructure = set(sum(substructure, []))
            flattened_substructures.append([x for (y, x)
                                           in sorted(zip(
                                               [self.levels[i] for i in
                                                flattened_substructure],
                                            flattened_substructure))])
        body_order = []

        # sort indices by level in the spanning tree
        for node in [x for (y, x) in sorted(zip([self.levels[i] for i in self.I], self.I))]:
            # create and expand short list simultaneously
            if not any([node in substructure for substructure in flattened_substructures]):
                body_order.append(node)
            else:
                for substructure in flattened_substructures:
                    if node is substructure[0]:
                        substructure_subset = [j for j in substructure if j not in body_order]
                        body_order = body_order + substructure_subset
        self.body_order = body_order
        self.joint_order = [self.parent_joints[i] for i in body_order][1:]


if __name__ == "__main__":
    # 0 1 2 3 4 5 6 7 8 9 10
    P = [-1, 0, 1, 2, 3, 4, 7, 4, 2, 2, 3,
         3, 3, 4, 4, 14, 8, 9, 12, 5, 15]
    loop_terminators = [[16, 10], [17, 11], [18, 13], [19, 6], [20, 6]]
    parent_joints = [None, 'A', 'B', 'C', 'D', 'E', 'T', 'F', 'H', 'I',
                     'P', 'Q', 'J', 'R', 'K', 'O', 'L', 'M', 'N', 'G', 'S']
    spanning_tree = SpanningTree(P, loop_terminators, parent_joints)
    print('P =', spanning_tree.P)
    print('L =', spanning_tree.L)
    print('R =', spanning_tree.R)
    print('levels =', spanning_tree.levels)
    print('substructures =')
    for index, substructure in enumerate(spanning_tree.substructures):
        print(index, substructure)
    print('revised body order =', spanning_tree.body_order)
    print('revised joint order =', spanning_tree.joint_order)
