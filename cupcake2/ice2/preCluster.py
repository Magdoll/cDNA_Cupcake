__author__ = 'etseng@pacb.com'

import random
from collections import defaultdict

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

    #def add_tucked(self, seqids):
    #    self.tucked += seqids


class preClusterSet:
    def __init__(self):
        self.S = {}
        self.seq_map = {}
        self.max_cid = 0
        self.tucked = defaultdict(lambda: []) # qID --> list of sIDs it is tucked under (can have multiple)


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
        seqid could already be in a cluster or tucked
        match could already be in a cluster or tucked
        so must handle this properly

        1. both in clusters --> combine
        2. one in cluster, other is not --> add other to the cluster
        3. neither in clusters --> create new cluster
        """
        if seqid in self.tucked:
            if match in self.tucked:
                # both is tucked, ignore
                return
            else:
                if match in self.seq_map:
                    # seqid is tucked, match is NOT tucked in 'c'
                    # make all members in 'c' tucked and removing it from S
                    cid = self.seq_map[match]
                    for x in self.S[cid].members:
                        self.tucked[x] = self.tucked[seqid]
                        del self.seq_map[x]
                    del self.S[cid]
                else:
                    self.tucked[match] = self.tucked[seqid]
                return
        else:
            if match in self.tucked:
                # seqid is not tucked, match is tucked
                if seqid in self.seq_map:
                    cid = self.seq_map[seqid]
                    for x in self.S[cid].members:
                        self.tucked[x] = self.tucked[match]
                        del self.seq_map[x]
                    del self.S[cid]
                else:
                    self.tucked[seqid] = self.tucked[match]
                return
            else:
                # both is untucked
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

    def add_seqid_contained(self, seqid, match):
        """
        seqid is contained in match
        1. if neither are in cluster, make "match" a cluster by itself and tuck seqid
        2. if seqid is in 'c' but not match, make "match" a cluster and move all in 'c' to tucked
        3. if seqid not in cluster but match in 'c', simply tuck seqid
        4. if both are in cluster, c1 and c2, move all in 'c1' to tucked

        NOTE: sometimes may happen that both are in same cluster! beware!
        """
        if seqid in self.tucked:
            if match in self.tucked:
                # both is tucked, ignore
                return
            elif match in self.seq_map:
                # seqid in tucked, match in 'c'
                # simply add to seqid tucked
                self.tucked[seqid].append(match)
            else:
                # seqid in tucked, match not in tucked or {S}
                # make matched a new cluster
                # tuck seqid to match
                self.add_new_cluster([match])
                self.tucked[seqid].append(match)
            return
        elif seqid in self.seq_map:
            # seqid in {S}
            cid1 = self.seq_map[seqid]
            if match in self.tucked:
                # seqid in 'c', match in tucked
                # tuck everything in same as where match is tucked
                for x in self.S[cid1].members:
                    self.tucked[x] += self.tucked[match]
                    del self.seq_map[x]
                del self.S[cid1]
            elif match in self.seq_map:
                # seqid and match both in clusters
                # need to handle cases in which both are in same cluster!
                if cid1 == self.seq_map[match]:
                    # ignore and not do anything if both in same cluster
                    return
                for x in self.S[cid1].members:
                    self.tucked[x].append(match)
                    del self.seq_map[x]
                del self.S[cid1]
            else:
                # seqid in 'c', match not in tucked or {S}
                # case (2), make match a cluster and move all in c1 to tucked
                self.add_new_cluster([match])
                for x in self.S[cid1].members:
                    self.tucked[x].append(match)
                    del self.seq_map[x]
                del self.S[cid1]
        else:
            # seqid not in tucked or {S}
            if match in self.tucked:
                # seqid not in anything, match in tucked
                self.tucked[seqid].append(match)
            elif match in self.seq_map:
                # seqid not in anything, match in 'c'
                self.tucked[seqid].append(match)
            else:
                # both are in nothing
                self.add_new_cluster([match])
                self.tucked[seqid].append(match)


    def reassign_members_by_cid(self, cid1, cid2):
        """
        Tuck all contents of cid1 into cid2
        (parallel to this, tucked should be updated in main program)
        """
        for seqid in self.S[cid1].members:
            self.seq_map[seqid] = cid2
        # now can safely delete cid1
        del self.S[cid1]



class preClusterSet2:
    def __init__(self):
        self.S = {}
        self.seq_map = {}
        self.max_cid = 0
        self.seq_stat = defaultdict(lambda: None) # seqid --> status ("None" is unmatched, "M" for matched, "T" for tucked)


    def add_new_cluster(self, members):
        self.S[self.max_cid] = preCluster(self.max_cid, members)
        for x in members:
            assert self.seq_stat[x] is None # SANITY: should not have been matched before
            self.seq_map[x] = self.max_cid
            self.seq_stat[x] = 'M'
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
        seqid could already be in a cluster or tucked
        match could already be in a cluster or tucked
        so must handle this properly

        1. both in clusters --> combine
        2. one in cluster, other is not --> add other to the cluster
        3. neither in clusters --> create new cluster
        """
        if self.seq_stat[seqid] == 'T':
            # seqid is tucked
            # case 1: match is tucked, do nothing
            # case 2: match is in 'c', tuck 'c'
            # case 3: match is unmatched, mark match as tucked
            if self.seq_stat[match] == 'M':
                self.tuck_cluster(self.seq_map[match])
            elif self.seq_stat[match] is None:
                self.seq_stat[match] = 'T'
        elif self.seq_stat[seqid] == 'M':
            # seqid is in 'c'
            # case 1: match is tucked, tuck 'c'
            # case 2: match is in 'd', if c!='d':combine
            # case 3: match is unmatched, add to 'c'
            cid1 = self.seq_map[seqid]
            if self.seq_stat[match] == 'T':
                self.tuck_cluster(cid1)
            elif self.seq_stat[match] == 'M':
                cid2 = self.seq_map[match]
                if cid1!=cid2:
                    self.combine_cluster(cid1, cid2)
            else:
                self.add_seqid_to_cluster_by_cid(match, cid1)
        else:
            # seqid is unmatched
            # case 1: match is tucked, also mark seqid as tucked
            # case 2: match is in 'c', add seqid to 'c'
            # case 3: match is unmatched, make a new cluster
            if self.seq_stat[match] == 'T':
                self.seq_stat[seqid] = 'T'
            elif self.seq_stat[match] == 'M':
                self.add_seqid_to_cluster_by_cid(seqid, self.seq_map[match])
            else:
                self.add_new_cluster([seqid, match])

    def add_seqid_to_cluster_by_cid(self, seqid, cid):
        assert seqid not in self.seq_map
        self.S[cid].add_member(seqid)
        self.seq_map[seqid] = cid
        self.seq_stat[seqid] = 'M'

    def add_seqid_to_cluster_by_match(self, seqid, match):
        self.add_seqid_to_cluster_by_cid(seqid, cid=self.seq_map[match])

    def tuck_cluster(self, cid):
        """
        Helper function: delete all members in 'cid', mark them as tucked ('T')
        """
        for x in self.S[cid].members:
            del self.seq_map[x]
            self.seq_stat[x] = 'T'
        del self.S[cid]

    def add_seqid_contained(self, seqid, match):
        """
        seqid is contained in match
        1. if neither are in cluster, make "match" a cluster by itself and tuck seqid
        2. if seqid is in 'c' but not match, make "match" a cluster and move all in 'c' to tucked
        3. if seqid not in cluster but match in 'c', simply tuck seqid
        4. if both are in cluster, c1 and c2, move all in 'c1' to tucked

        NOTE: sometimes may happen that both are in same cluster! beware!
        """
        if self.seq_stat[seqid] == 'T':
            # if seqid is tucked,
            # case 1: match is tucked, do nothing
            # case 2: match is in 'c', do nothing
            # case 3: match is unmatched, make it a cluster
            if self.seq_stat[match] is None:
                self.add_new_cluster([match])
        elif self.seq_stat[seqid] == 'M':
            # if seqid is in 'c',
            # case 1: match is tucked, mark everything 'c' as tucked, delete 'c'
            # case 2: match is in 'd', if c==d, ignore (weird case!); if c!=d, delete 'c'
            # case 3: match is unmatched, make it a cluster and delete 'c'
            cid1 = self.seq_map[seqid]
            if self.seq_stat[match] == 'T':
                self.tuck_cluster(cid1)
            elif self.seq_stat[match] == 'M':
                cid2 = self.seq_map[match]
                if cid1!=cid2:
                    self.tuck_cluster(cid1)
            else:
                self.add_new_cluster([match])
                self.tuck_cluster(cid1)
        else:
            # seqid is unmatched
            # case 1: match is tucked, mark seqid as tucked
            # case 2: match is in 'c', mark seqid as tucked
            # case 3: match is unmatched, make it a cluster, mark seqid as tucked
            self.seq_stat[seqid] = 'T'
            if self.seq_stat[match] is None:
                self.add_new_cluster([match])
