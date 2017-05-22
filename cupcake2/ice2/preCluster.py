__author__ = 'etseng@pacb.com'

import random

class preCluster:
    def __init__(self, cid, members=[]):
        self.cid = cid
        self.size = len(members)
        self.members = members
        self.tucked = []

    def add_member(self, seqid):
        self.size += 1
        self.members.append(seqid)

    def add_members(self, seqids):
        self.size += len(seqids)
        self.members += seqids

    def add_tucked(self, seqids):
        self.tucked += seqids


class preClusterSet:
    def __init__(self):
        self.S = {}
        self.seq_map = {}
        self.max_cid = 0


    def add_new_cluster(self, members):
        self.S[self.max_cid] = preCluster(self.max_cid, members)
        for x in members:
            self.seq_map[x] = self.max_cid
        self.max_cid += 1
        return self.max_cid-1

    def combine_cluster(self, cid1, cid2):
        """
        just make all members of cid2 into members of cid1
        """
        assert cid1!=cid2
        # update members for cid1
        self.S[cid1].add_members(self.S[cid2].members)
        # update seq_map
        for x in self.S[cid2].members:
            self.seq_map[x] = cid1
        # now can safely delete cid2
        del self.S[cid2]

    def add_seqid_match(self, seqid, match):
        """
        seqid could already be in a cluster
        match could already be in a cluster
        so must handle this properly

        1. both in clusters --> combine
        2. one in cluster, other is not --> add other to the cluster
        3. neither in clusters --> create new cluster
        """
        if seqid in self.seq_map:
            cid1 = self.seq_map[seqid]
            if match in self.seq_map:
                # both are in clusters, combine!
                cid2 = self.seq_map[match]
                if cid1 != cid2:
                    self.combine_cluster(cid1, cid2)
            else:
                self.add_seqid_to_cluster_by_cid(match, cid1)
        else:
            if match in self.seq_map:
                cid2 = self.seq_map[match]
                self.add_seqid_to_cluster_by_cid(seqid, cid2)
            else:
                # neither are in clusters, create a new one
                self.add_new_cluster([seqid, match])

    def add_seqid_to_cluster_by_cid(self, seqid, cid):
        assert seqid not in self.seq_map
        self.S[cid].add_member(seqid)
        self.seq_map[seqid] = cid

    def add_seqid_to_cluster_by_match(self, seqid, match):
        cid = self.seq_map[match]
        self.add_seqid_to_cluster_by_cid(seqid, cid)

    def reassign_members_by_cid(self, cid1, cid2):
        """
        Tuck all contents of cid1 into cid2
        (parallel to this, tucked should be updated in main program)
        """
        for seqid in self.S[cid1].members:
            self.seq_map[seqid] = cid2
        # now can safely delete cid1
        del self.S[cid1]