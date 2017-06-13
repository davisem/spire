import tensorflow as tf
import numpy as np
import matplotlib.pyplot as plt


num_hashes = 400
num_buckets = 2**32
cpu = "/cpu:0"
window_length = 1000
kmer_length = 8
prob_error = 0.4

minhash = range(num_hashes)

with tf.device(cpu):
    kmers = tf.placeholder( tf.string, shape=[None], name="kmers")
    original_hashed_kmers = tf.string_to_hash_bucket_fast(kmers, num_buckets, name="hashed_kmers")
    hashed_kmers = original_hashed_kmers
    for i in range(num_hashes-1):
        minhash[i] = tf.argmin( hashed_kmers )
        hashed_kmers = tf.scalar_mul( 0x1234567887654321, hashed_kmers)
        hashed_kmers = 0x1234567887654321 + hashed_kmers
    minhash[num_hashes-1] = tf.argmin( hashed_kmers )
 
def MinHash(A):
    with tf.Session(config=tf.ConfigProto(log_device_placement=False)) as sess:
        my_minhash, A_strings = sess.run( [minhash,kmers],  feed_dict={kmers: A} )
    return A_strings[my_minhash]

def Jaccard(A,B):
    hA = MinHash(A.keys())
    hB = MinHash(B.keys())
    return sum(hA==hB)*1.0/num_hashes

def Jaccardt(A,B,t):
    s = 0
    hA = MinHash(A.keys())
    hB = MinHash(B.keys())
    for i in range(num_hashes) :
        if (hA[i] == hB[i]) and ( abs(A[hA[i]] - B[hB[i]]) < t ):
            s += 1
    return s*1.0/num_hashes     


def JaccardRange( A, B, tmin, tmax ):
    jac = range(tmax) 
    hA = MinHash(A.keys())
    hB = MinHash(B.keys())
    for t in range(tmin,tmax) : 
        s = 0
        for i in range(num_hashes) :
            if (hA[i] == hB[i]) and ( abs(A[hA[i]] - B[hB[i]]) < t ):
                s += 1
        jac[t] = s*1.0/num_hashes
        print t,s*1.0/num_hashes
    return jac

def CreateWindow( wl ):
    window = ""
    r = np.random.randint(4, size=wl)
    for i in range(wl):
        if r[i] == 0:
            window += 'G'
        elif r[i] == 1:
            window += 'T'
        elif r[i] == 2:
            window += 'C'
        else:
            window += 'A' 
    return window


def PerturbWindow ( original_window, prob ):
    window = ""
    r = np.random.rand(len(original_window))
    for i in range ( len(original_window)):
        if r[i] < prob :
            rr = np.random.randint(6)
            if rr == 0:
                window += 'G'
            elif rr == 1:
                window += 'T'
            elif rr == 2:
                window += 'C'
            elif rr == 3:
                window += 'A'
            elif rr == 4:
                rrr = np.random.randint(4)
                if rrr == 0:
                    window += 'G'
                elif rrr == 1:
                    window += 'T'
                elif rrr == 2:
                    window += 'C'
                elif rrr == 3:
                    window += 'A'
                i -= 1
        else:
            window += original_window[i]
    return window           

def CreateKMERPairs ( w, k ):
    kmerpairs = {}
    for i in range ( len(w) - k):
        for j in range (i+1,len(w) - k ):
            a = w[i:i+k]
            b = w[j:j+k]
            d = j - i
            ab= a+","+b+":"
            if ab in kmerpairs:
                kmerpairs[ab] = min(d,kmerpairs[ab])
            else :
                kmerpairs[ab] = d
    return kmerpairs


def CreateKMERs ( w, k ):
    kmerpairs = {}
    for i in range ( len(w) - k):
            a = w[i:i+k]
            ab= a+":"
            kmerpairs[ab] = 0
    return kmerpairs




def CreateKMERPairsThreshold ( kmers, t ):
    newkmers = {}
    for k in kmers.keys():
        if kmers[k] < t:
            newkmers[k+"0"] = 0
        else :
            newkmers[k+"1"] = 1
    return newkmers




w1 = CreateWindow( window_length )
kmerpairs1 = CreateKMERPairs ( w1, kmer_length )
kmers1 = CreateKMERs ( w1, 2*kmer_length )
w2 = PerturbWindow( w1, prob_error )
kmerpairs2 = CreateKMERPairs ( w2, kmer_length )
kmers2 = CreateKMERs ( w2, 2*kmer_length )
w3 = CreateWindow( window_length )
kmerpairs3 = CreateKMERPairs ( w3, kmer_length )
kmers3 = CreateKMERs ( w3, 2*kmer_length )

jacpair = JaccardRange(kmerpairs1, kmerpairs2, 0, window_length-kmer_length )
jac = Jaccard(kmers1, kmers2 )
jacpairrandom = JaccardRange(kmerpairs1, kmerpairs3, 0, window_length-kmer_length )
jacrandom = Jaccard(kmers1, kmers3 )


f = plt.figure()
plt.plot(jacpair, 'r', label="kmer pairs error="+str(prob_error))
plt.plot(jacpairrandom, 'g', label="kmer pairs, random")
plt.plot([0,window_length],[jac,jac],'b', label="2kmer")
plt.plot([0,window_length],[jacrandom,jacrandom],'y',label="2kmer, random")
plt.legend()
plt.xlabel('threshold', fontsize=14)
plt.ylabel('IoU, Jaccard', fontsize=14)
plt.title('Intersection over Union for error='+str(prob_error)+' ,k= '+str(kmer_length))
plt.show()
f.savefig("graph"+str(prob_error)+"."+str(kmer_length)+".pdf", bbox_inches='tight')
