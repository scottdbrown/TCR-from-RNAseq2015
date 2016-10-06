import sys
import os
import random

'''
TCR read simulator
@author: sbrown

Generates a fastq file of alpha and beta TCRs
Input: location of TCR gene references, output file, number of TCRs for each chain

'''

ref = sys.argv[1]
out = sys.argv[2]
NUM_TCR = int(sys.argv[3])





## find files in `ref`:

av = os.path.join(ref, "TRAV.fa")
aj = os.path.join(ref, "TRAJ.fa")
ac = os.path.join(ref, "TRAC.fa")

bv = os.path.join(ref, "TRBV.fa")
bd = os.path.join(ref, "TRBD.fa")
bj = os.path.join(ref, "TRBJ.fa")
bc = os.path.join(ref, "TRBC.fa")

''' CDR3 insertion/deletion probabilities '''
''' from Warren et al. Profiling model T-cell metagenomes with short reads. Table 1 '''
## pVDel is probability of number of bases of 3' V deletion
## [number of bases, probability]
pVDel = [[0,0.194],[1,0.160],[2,0.098],[3,0.118],[4,0.160],[5,0.118],[6,0.070],[7,0.045],[8,0.022],[9,0.008],[10,0.006]]
pJDel = [[0,0.209],[1,0.123],[2,0.122],[3,0.104],[4,0.117],[5,0.123],[6,0.086],[7,0.056],[8,0.031],[9,0.019],[10,0.010]]
pCDR3Add = [[1,0.006],[3,0.006],[4,0.029],[5,0.017],[6,0.017],[7,0.052],[8,0.063],[9,0.069],[10,0.080],[11,0.057],[12,0.075],[13,0.080],[14,0.086],[15,0.098],[16,0.046],[17,0.052],[18,0.029],[19,0.034],[20,0.023],[21,0.023],[22,0.011],[23,0.006],[25,0.011],[26,0.011],[27,0.006],[28,0.011]]

''' from Saada et al. Models for antigen receptor gene rearrangment: CDR3 length. '''
## frequency of each base as non-templated addition.
bases = [["A",0.20],["T",0.20],["C",0.22],["G",0.38]]

''' Functions '''
def createRefDict(file):
    '''Returns dictionary of sequences'''
    dic = {}
    gene = ""
    id = 0
    for line in open(file, "r"):
        if line[:1] == ">":
            gene = line[1:].rstrip()
            id += 1
        else:
            dic[id] = [gene, line.rstrip()]
    return dic


def weighted_choice(items):
    '''items is a list of tuples in the form (item, weight)'''
    weight_total = sum((item[1] for item in items))
    n = random.uniform(0, weight_total)
    for item, weight in items:
        if n < weight:
            return item
        n = n - weight
    return item


def containStop(seq):
    ''' recursively checks for presence of in-frame stop codon'''
    stops = ["TAA","TAG","TGA"]
    ## base case:
    if len(seq) < 3:
        return False
    ## check for stop codons
    elif seq[:3] in stops:
        return True
    else:
        return containStop(seq[3:])




''' END Functions '''


## Create dictionaries:

ava = createRefDict(av)
aja = createRefDict(aj)
aca = createRefDict(ac)

bva = createRefDict(bv)
bda = createRefDict(bd)
bja = createRefDict(bj)
bca = createRefDict(bc)


outfile = open(out, "w")

random.seed(1)  ## define randomness for now

for i in range(1,NUM_TCR+1):

    inFrame = False
    orf = False

    while not (inFrame and orf):

        seq = ""
        vDel = 0
        jDel = 0
        inFrame = False
        orf = False

        ## make a TRA

        ## get random index for each.
        rav = random.randint(1,len(ava))
        raj = random.randint(1,len(aja))
        rac = random.randint(1,len(aca))

        ## get gene sequences
        vGene = ava[rav][1]
        jGene = aja[raj][1]
        cGene = aca[rac][1]

        ## determine number of bases to delete and add (using weighted scores)
        vDel = weighted_choice(pVDel)
        jDel = weighted_choice(pJDel)
        cdr3Add = weighted_choice(pCDR3Add)

        ## determine if will result in in-frame. if not, try again.
        if (len(vGene) - vDel + cdr3Add + len(jGene) - jDel + len(cGene)) % 3 == 0:
            inFrame = True
        else:
            inFrame = False
            continue    ## no sense in doing mutations...try again.

        ## do somatic mutations.
        ## Delete bases off 3' end of V:
        if vDel != 0:
            vGene = vGene[:-vDel]
        ## Delete bases off 5' end of J:
        if jDel != 0:
            jGene = jGene[jDel:]
        ## Add CDR3
        addedBases = ""
        for j in range(1,cdr3Add+1):
            ## Add bases between V and J
            addedBases += weighted_choice(bases)
        ## check for stop codons
        ## join all pieces into a spliced mRNA:
        seq = vGene + addedBases + jGene + cGene

        orf = not containStop(seq[:-3])     #cut off last codon as it should be a stop.


    id = "TRA" + str(i) + ":" + ava[rav][0] + ":" + str(vDel) + ":" + addedBases + ":" + str(jDel) + ":" + aja[raj][0] + ":" + aca[rac][0] + ":" + str(len(seq))

    outfile.write("@" + id + "\n")
    outfile.write(seq + "\n")
    outfile.write("+\n")
    outfile.write("J"*(len(seq)))
    outfile.write("\n")


    ## Make T cell receptor beta

    inFrame = False
    orf = False

    while not (inFrame and orf):
        seq = ""
        vDel = 0
        jDel = 0
        inFrame = False
        orf = False

        ## make a TRB
        #print "DEBUG: len(bva)=" + str(len(bva)) + "\n"
        rbv = random.randint(1,len(bva))
        rbd = random.randint(1,len(bda))
        rbj = random.randint(1,len(bja))
        rbc = random.randint(1,len(bca))

        vGene = bva[rbv][1]
        dGene = bda[rbd][1]
        jGene = bja[rbj][1]
        cGene = bca[rbc][1]

        ## determine number of bases to delete and add (using weighted scores)
        vdDel = weighted_choice(pVDel)
        djDel = weighted_choice(pJDel)

        # split the vdDel between 3' of V and 5' of D
        vDel = random.randint(0,vdDel)
        d5Del = vdDel - vDel

        # split the djDel between 3' of D and 5' of J
        jDel = random.randint(0,djDel)
        d3Del = djDel - jDel

        cdr3Add = weighted_choice(pCDR3Add)
        vdAdd = random.randint(0,cdr3Add)
        djAdd = cdr3Add - vdAdd


        ## determine if will result in in-frame. if not, try again.
        if (len(vGene) - vdDel + cdr3Add + +len(dGene) + len(jGene) - djDel + len(cGene)) % 3 == 0:
            inFrame = True
        else:
            inFrame = False
            continue    ## no sense in doing mutations...try again.

        ## do somatic mutations.
        ## Delete bases off 3' end of V:
        if vDel != 0:
            vGene = vGene[:-vDel]
        ## Delete bases off 5' end of D:
        if d5Del != 0:
            dGene = dGene[d5Del:]
        ## Delete bases off 3' end of D:
        if d3Del != 0:
            dGene = dGene[:-d3Del]
        ## Delete bases off 5' end of J:
        if jDel != 0:
            jGene = jGene[jDel:]
        ## Add CDR3
        vdAddedBases = ""
        for j in range(1,vdAdd+1):
            ## Add bases between V and D
            vdAddedBases += weighted_choice(bases)
        djAddedBases = ""
        for j in range(1,djAdd+1):
            ## Add bases between D and J
            djAddedBases += weighted_choice(bases)
        ## check for stop codons
        ## join all pieces into a spliced mRNA:
        seq = vGene + vdAddedBases + dGene + djAddedBases + jGene + cGene

        orf = not containStop(seq[:-3])     #cut off last codon

    #print "DEBUG: rbv: " + str(rbv) + "\n"

    id = "TRB" + str(i) + ":" + bva[rbv][0] + ":" + str(vDel) + ":" + vdAddedBases + ":" + str(d5Del) + ":" + bda[rbd][0] + ":" + str(d3Del) + ":" + djAddedBases + ":" + str(jDel) + ":" + bja[rbj][0] + ":" + bca[rbc][0] + ":" + str(len(seq))

    outfile.write("@" + id + "\n")
    outfile.write(seq + "\n")
    outfile.write("+\n")
    outfile.write("J"*(len(seq)))
    outfile.write("\n")


outfile.close()
print "done."
