import tensorflow as tf
import numpy as np
import matplotlib.pyplot as plt


num_hashes = 400
num_buckets = 2**32
cpu = "/cpu:0"
window_length = 1000
kmer_length = 8
prob_error = 0.2
num_trials = 100


minhash = range(num_hashes)

with tf.device(cpu):
    r = np.random.randint(num_buckets,2*num_buckets,[num_hashes,2])
    kmers = tf.placeholder( tf.string, shape=[None], name="kmers")
    for i in range(num_hashes):
        minhash[i] = tf.argmin( tf.string_to_hash_bucket_strong(kmers, num_buckets, [r[i,0],r[i,1]]) )
 
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


jacpair = []
jac = []
jacpairrandom = []
jacrandom = []
for i in range(num_trials):
    print i
    w1 = CreateWindow( window_length )
    kmerpairs1 = CreateKMERPairs ( w1, kmer_length )
    kmers1 = CreateKMERs ( w1, 2*kmer_length )
    w2 = PerturbWindow( w1, prob_error )
    kmerpairs2 = CreateKMERPairs ( w2, kmer_length )
    kmers2 = CreateKMERs ( w2, 2*kmer_length )
    w3 = CreateWindow( window_length )
    kmerpairs3 = CreateKMERPairs ( w3, kmer_length )
    kmers3 = CreateKMERs ( w3, 2*kmer_length )

    jacpair.append(JaccardRange(kmerpairs1, kmerpairs2, 0, window_length-kmer_length ))
    jac.append(Jaccard(kmers1, kmers2 ))
    jacpairrandom.append(JaccardRange(kmerpairs1, kmerpairs3, 0, window_length-kmer_length ))
    jacrandom.append(Jaccard(kmers1, kmers3 ))

jacpair_mean = np.mean(jacpair,0)
jacpair_std = np.std(jacpair,0)
jac_mean = np.mean(jac)
jac_std = np.std(jac)
jacpairrandom_mean = np.mean(jacpairrandom,0)
jacpairrandom_std = np.std(jacpairrandom,0)
jacrandom_mean = np.mean(jacrandom)
jacrandom_std = np.std(jacrandom)


f = plt.figure()
plt.plot(jacpair_mean, 'r', label="kmer pairs, error="+str(prob_error))
plt.plot(jacpair_mean+jacpair_std, 'r')
plt.plot(jacpair_mean-jacpair_std, 'r')
plt.plot(jacpairrandom_mean, 'g', label="kmer pairs, random")
plt.plot(jacpairrandom_mean+jacpairrandom_std, 'g')
plt.plot(jacpairrandom_mean-jacpairrandom_std, 'g')
plt.plot([0,window_length],[jac_mean,jac_mean],'b', label="2kmeri, error="+str(prob_error))
plt.plot([0,window_length],[jac_mean+jac_std,jac_mean+jac_std],'b')
plt.plot([0,window_length],[jac_mean-jac_std,jac_mean-jac_std],'b')
plt.plot([0,window_length],[jacrandom_mean,jacrandom_mean],'y',label="2kmer, random")
plt.plot([0,window_length],[jacrandom_mean+jacrandom_std,jacrandom_mean+jacrandom_std],'y')
plt.plot([0,window_length],[jacrandom_mean-jacrandom_std,jacrandom_mean-jacrandom_std],'y')
plt.legend()
plt.xlabel('threshold', fontsize=14)
plt.ylabel('IoU, Jaccard', fontsize=14)
plt.title('Intersection over Union for error='+str(prob_error)+' ,k= '+str(kmer_length))
plt.show()
f.savefig("plots/graph"+str(prob_error)+"."+str(kmer_length)+".pdf", bbox_inches='tight')
