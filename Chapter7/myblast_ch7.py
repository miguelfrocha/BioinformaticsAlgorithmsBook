
def read_database (filename):
    f = open (filename)
    db = []
    for line in f:
        db.append(line.rstrip())
    f.close()
    return db

def build_map (query, w):
    res = {}
    for i in range(len(query)-w+1):
        subseq = query[i:i+w]
        if subseq in res:
            res[subseq].append(i)
        else:
            res[subseq] = [i]
    return res

def get_hits (seq, m, w):
    res = []  # list of tuples
    for i in range(len(seq)-w+1):
        subseq = seq[i:i+w]
        if subseq in m:
            l = m[subseq]
            for ind in l:
                res.append( (ind,i) )
    return res
    
def extends_hit (seq, hit, query, w):
    stq, sts = hit[0], hit[1]
    ## move forward
    matfw = 0       
    k=0
    bestk = 0
    while 2*matfw >= k and stq+w+k < len(query) and sts+w+k < len(seq):
        if query[stq+w+k] == seq[sts+w+k]: 
            matfw+=1
            bestk = k+1
        k += 1
    size = w + bestk
    ## move backwards
    k = 0
    matbw = 0   
    bestk = 0
    while 2*matbw >= k and stq > k and sts > k:
        if query[stq-k-1] == seq[sts-k-1]: 
            matbw+=1
            bestk = k+1
        k+=1       
    size += bestk
    
    return (stq-bestk, sts-bestk, size, w+matfw+matbw)


def hit_best_score(seq, query, m, w):
    hits = get_hits(seq, m, w)
    bestScore = -1.0
    best = ()
    for h in hits:
        ext = extends_hit(seq, h, query, w)
        score = ext[3]
        if score > bestScore or (score== bestScore and ext[2] < best[2]):
            bestScore = score
            best = ext
    return best
 
 
def best_alignment (db, query, w):
    m = build_map(query, w)
    bestScore = -1.0
    res = (0,0,0,0,0)
    for k in range(0,len(db)):
        bestSeq = hit_best_score(db[k], query, m, w)
        if bestSeq != ():
            score = bestSeq[3]  
            if score > bestScore or (score== bestScore and bestSeq[2] < res[2]):
                bestScore = score
                res = bestSeq[0], bestSeq[1], bestSeq[2], bestSeq[3], k
    if bestScore < 0: return ()
    else: return res


 
# test

db = read_database("seqBlast.txt")
query = "gacgcctcgcgctcgcgcgctgaggcaaaaaaaaaaaaaaaaaaaatcggatagctagctgagcgctcgatagcgcgttcgctgcatcgcgtatagcgctgaagctcccggcgagctgtctgtaaatcggatctcatctcgctctatcct"
r = best_alignment(db, query, 11)
print(r)

query2 = "cgacgacgacgacgaatgatg"
r = best_alignment(db, query2, 11)
print(r)