#Naam: Sil van Appeldoorn
#Datum: 20-01-2018
#Deze afvinkopdracht is gemaakt met de hulp van Bram Bosch, omdat ik net deze weken had gemist door persoonlijke omstandigheden en ik deze afvink heel moeilijk vond heeft hij mij geholpen met deze afvink.
#We hebben hem samen stap voor stap uitgewerkt zodat ik het snapte en wist wat ik deed.

import re
import time

def main():
    bestand_naam = "coli_dna.fasta"
    
    bestand = open(bestand_naam)
    
    start_time = time.time()
    uitzonderingIT = is_dna(bestand)
    print("--- %s seconds door te ittereren ---" %  (time.time() - start_time))

    bestand.close()
    bestand = open(bestand_naam)
    
    start_time = time.time()
    uitzonderingRE = check_dna(bestand)
    print("--- %s seconds met regular expresions ---" % (time.time() -  start_time))

    print("Niet-DNA regels met regular expresion regels",uitzonderingRE)
    print("Niet-DNA regels door te ittereren", uitzonderingIT)

    bestand.close()
    bestand = open(bestand_naam) 

    
    if uitzonderingRE == 0 and uitzonderingIT == 0:
        bestandstype = "dna"
    else:
        bestandstype = "eiwit"
    
    if bestandstype == "eiwit":

        seqs = lees_eiwit(bestand)
        begin = []
        plek_consensus = []
        for header in seqs:
            begin, plek_consensus = find_p53(plek_consensus, begin, seqs, header)
        i=0
        
        for waarde in plek_consensus:
            print("De header waarin het consensus patroon staat is:",plek_consensus[i])
            print("In deze header staat het consensus patroon op positie:",begin[i])
            i += 1
            print("\n")
        for header in seqs:
            find_aa(header)
    else:
        
        sequentie = lees_dna(bestand)
        if (re.search("TTGACA[ATGC]{10,50}TATA[ATCG]{10}(ATG)([ATGC]{3})*(TGA|TAA|TAG)",sequentie)) == None:
            print("De promotor staat niet in de sequentie")
        else:
            print("TTGACA[ATGC]{10,50}TATA[ATCG]{10}(ATG)([ATGC]{3})*(TGA|TAA|TAG)",sequentie)
         
def lees_dna(bestand):
    sequentie = ""
    for regel in bestand:
        if regel.startswith(">"):
            pass
        else:
            sequentie += regel        
    return sequentie
    
    
def check_dna(bestand):

    uitzonderingRE = 0
    for regel in bestand:
        regel = regel.strip()
        dna = re.match("^[ATCG]*$",regel)
        if dna == None:
            if regel.startswith(">"):
                pass
            else:
                uitzonderingRE +=1           
    return uitzonderingRE
            
    

def lees_eiwit(bestand):    
    headers = "a"
    seqs = {}
    string = ""
    
    
    for line in bestand:
        line = line.strip()
        if line.startswith(">"):
            if len(headers) == 1:
                line = line.split()
                headers = line[0]
                string = ""
                
            else:
                line = line.split(" ")
                seqs[headers] = string
                headers = line[0]
                string = ""
                
        else:
            string += line
    seqs[headers] = string

    return seqs
    
def find_p53(plek_consensus, begin, seqs, header):
    
    consensus = re.search("MCNSSC[MV]GGMNRR",seqs[header])
    if consensus != None:
        begin.append(consensus.start())
        plek_consensus.append(header)
    return begin, plek_consensus
    
    
def find_aa(eiwit):
    if re.match("^[ARNDCFQEGHILKMPSTWYV]*$",eiwit):
        print("De sequentie bevat non-aminozuur tekens")

def is_dna(bestand):

    dna = []
    a=0
    t=0
    c=0
    g=0
    uitzonderingIT = 0
    for regel in bestand:
        regel = regel.strip()
        if regel.startswith(">"):
            pass
        else: 
            a = regel.count("A")
            t = regel.count("T")    
            c = regel.count("C")
            g = regel.count("G")
            if (a + t + c + g) == len(regel):
                dna.append(True)
            else:
                dna.append(False)
                
    
    for item in dna:
        if item == False:
            uitzonderingIT += 1
        
    return uitzonderingIT

main()
