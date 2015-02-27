 #       __  __  ___  _   _  ___  ____  
 #      |  \/  |/ _ \| \ | |/ _ \|  _ \ 
 #      | |\/| | | | |  \| | | | | | | |
 #      | |  | | |_| | |\  | |_| | |_| |
 #      |_|  |_|\___/|_| \_|\___/|____/ 
 #      robert harry nicodemus williams
 #
 #      Mayday 2009 -- 
 #
 #   GOALS:
 #    (Trivial)   [ ] 2 + 2 = 4
 #    (Easy)      [ ] Hello, World!
 # (Pretty easy?) [ ] Fibonacci
 #    (hard?)     [ ] Brainfuck


import random, copy, string, sys

CODONLENGTH = 3

class Code:
  def __init__(self):
    global basealphabet
    basealphabet = ("A","T","G","C")
    global operatorlist
    operatorlist = ("NULL","FOOO","LYSE","META","CATA","PHOS","DEPH","COPH","BIND","BEND","BACK","TRAN","BEGN","STOP","SENS","JESU")

    # reminders for what these are tentatively supposed to do:
    #
    # NULL: 0
    # FOOO: a variant of NULL
    #
    # LYSE: Destroy the protein /(it is bound to)
    #
    # META + number: Increment chemical_{number} by phosphorylation status
    # CATA + number: Decrement blah blah
    #
    # PHOS: increment phosphorylation of protein /(it is bound to)
    # DEPH: decrement blah blah
    # COPH: COPH n -> phosphorylate if register n is nonzero
    #
    # BIND: mark start of binding sequence
    # BEND: mark end of binding sequence
    # BACK: unbind proteins and attempt to bind to first protein with binding sequence following.
    #
    # TRAN: this one starts transcription of proteins - protein under construction is automatically BACK-bound to creating protein and should prob. be BACK-unbound after 
    # BEGN: a new protein begins here: make operons (not yet) with BIND *brbrbrbrbr BEND BEGN PROT1 PROT1 BEGN PROT2 PROT2 STOP
    # STOP: transcription ceases entirely
    #
    # SENS: "SENS n" : register[n] -> phosphorylation
    # JESU: "JESU n" : n -> register[phosphorylation]
    #--------------------------------------------------------------------
    #proposition 
    #TAKE: TAKE n would halt computation and ask for a value to put into memory address n, then continue. or maybe phos status
    #KINA: ph(-1) self, ph(+1) bound, if unbound, ph(0)
    #
    #
    # Kadmon protein: transcribes *tatgag protein into existence, releases protein (BACK) then kills itself
    # TRAN BIND *tatgag BEND BACK LYSE
    #
    # genome then is like:
    # NULL NULL NULL NULL NULL BIND *tatgag BEND BEGN FOOO FOOO FOOO STOP NULL NULL NULL

    #this is painfully pre v1.0,
    #           ("NULL","FOOO","LYSE","META","CATA","PHOS","DEPH","COPH","BIND","BEND","BACK","TRAN","BEGN","STOP","SENS","JESU")
   #canonical = ["AAA" ,"TTT" ,"CAG" ,"TAC" ,"CAT" ,"CCC" ,"GGG" ,"GAC", "TAA" ,"ATA" ,"AAT" ,"CGG" ,"GCG" ,"GGC", "ACG", "ACT" ]
    
    #for turning ???? back into DNA
    global ballscodons
    ballscodons = ("ACA","ACC")
    
    global codebook
    #pre v1.0 again
    #codebook = dict(zip(canonical,operatorlist))
    
    #operatorlist = ("NULL","FOOO","LYSE","META","CATA","PHOS","DEPH","COPH","BIND","BEND","BACK","TRAN","BEGN","STOP","SENS","JESU")
    
    codebook = {}

    codebook["AAA"] = "NULL";codebook["AAC"] = "NULL";codebook["AAG"] = "NULL"
    
    codebook["TTA"] = "FOOO";codebook["TTT"] = "FOOO";codebook["TTC"] = "FOOO";codebook["TTG"] = "FOOO"

    codebook["CAC"] = "LYSE";codebook["CAG"] = "LYSE"

    codebook["TAG"] = "META";codebook["TAC"] = "META"

    codebook["CAT"] = "CATA";codebook["CAA"] = "CATA"

    codebook["CCA"] = "PHOS";codebook["CCT"] = "PHOS";codebook["CCC"] = "PHOS";codebook["CCG"] = "PHOS"

    codebook["GGG"] = "DEPH"

    codebook["GAA"] = "COPH";codebook["GAT"] = "COPH";codebook["GAC"] = "COPH";codebook["GAG"] = "COPH"

    codebook["TAA"] = "BIND";codebook["TAT"] = "BIND" 
                         
    codebook["ATT"] = "BEND";codebook["ATC"] = "BEND";codebook["ATA"] = "BEND"


    codebook["AAT"] = "BACK"

    codebook["CGA"] = "TRAN";codebook["CGT"] = "TRAN";codebook["CGC"] = "TRAN";codebook["CGG"] = "TRAN"

    codebook["GCA"] = "BEGN";codebook["GCT"] = "BEGN";codebook["GCG"] = "BEGN";codebook["GCC"] = "BEGN"

    codebook["GGC"] = "STOP";codebook["GGA"] = "STOP";codebook["GGT"] = "STOP"

    codebook["ACG"] = "SENS"

    codebook["ACT"] = "JESU"

    # for errors
    default = "default"
    codebook["default"] = "????"

  def codon_to_integer(self,codon):
    """Converts a given codon into an unsigned integer.
       Each member of the alphabet is given a number, which gives each
       codon an n-ary value, which can be trivially converted to decimal"""
    #for 0-length inputs
    if len(codon) == 0:
      return 0
    
    numbase = len(basealphabet)
    indexdict = dict(zip(basealphabet,range(numbase)))
    
    #print codon
    numberlist = []
    integer = 0
    for letter in codon:
      numberlist.append(str(indexdict[letter]))
    #print numberlist
    integer = int("".join(numberlist),numbase)
    return integer
    
  def integer_to_codon(self,integer):
   """Takes an integer and converts it to a codon of approved length.
      If given an integer greater than is possible to express in the 
      base used, this function returns 0. If given an integer too low
      it will pad it with whatever is mapped to 0. """
   integer = int(integer)
   numbase = len(basealphabet)
   
   maxint = 0
   for i in range(numbase):
     maxint = maxint + (numbase**(numbase-i))
   if integer > maxint:
     return 0
   indexdict = dict(zip(range(numbase),basealphabet))
   
   wantednum = []
   for i in range(numbase-1):
     wantednum.append(indexdict[integer%numbase])
     integer = integer/numbase
   wantednum.reverse()
   wantedcodon = "".join(wantednum)
   return wantedcodon
   
  def codon_to_operator(self,codon):
    """Takes a codon, returns its operator"""
    operator = ""
    operator = codebook.get(codon)
    if operator == None:
      operator = "????"
    return operator

  def operator_to_codon(self,operator):
    """Takes an operator, returns the codon"""
    codon = ""
    for k in codebook:
      if codebook[k] == operator:
        return k
    raise ValueError

  def tokenize_dna(self,sequence):
    """Takes a list of strings in the alphabet and returns a tokenised version of it, i.e:
       CODONLENGTH length codons, with noncoding sequences put into large tokens, also dealing
       with fragments and stop codons"""
    pointer = 0
    tokens = []
    oncodons = ("ATA","ATT","ATC")  #(self.operator_to_codon("BEND"))
    offcodons = ("TAA","TAT")  #(self.operator_to_codon("BIND"))
    
    while pointer < len(sequence):
      ## if there's more than CODONLENGTH bases left in the sequence...
      if len(sequence) - pointer >= CODONLENGTH:

        subpointer = 0
        tempcodon = []
        while len(tempcodon) < CODONLENGTH:
          tempcodon.extend(sequence[pointer+subpointer])
          subpointer = subpointer + 1

        newcodon = "".join(tempcodon)
        
        # loop for binding sequences
        if newcodon in offcodons:
          temp = []
          tokens.append(newcodon)
          # until we reach a debinding sequence
          while newcodon not in oncodons:

            subpointer = 0
            tempcodon = []

            while len(tempcodon) < CODONLENGTH:
            # has to do with BIND with no closing BEND i think
              tempcodon.extend(sequence[pointer+subpointer])
              subpointer = subpointer + 1 
              
            newcodon = "".join(tempcodon)

            temp.append(newcodon)
            pointer = pointer + CODONLENGTH
          # when we do, we can finally dump all the codons we kept into a token, removing the signal
          # codons, of course!
          temp = temp[1:-1]
          temp = "".join(temp)
          # a hunch - okay this works, I don't know why. The tokeniser ignores the codon directly after BEND
          # unless you move this back a codon - presumably it goes one too far at some final step
          pointer = pointer - CODONLENGTH
          tokens.append(temp)
      
      # on the other hand, if there's less than three left, there's a fragment to deal with later...
      if len(sequence) - pointer < CODONLENGTH:
       
       fraglength = (len(sequence) - pointer)
       subpointer = 0
       tempcodon = []
       while len(tempcodon) <= fraglength:
         tempcodon.extend(sequence[pointer+subpointer])
         subpointer = subpointer + 1
       
       newcodon = "".join(tempcodon)
       
       tokens.append(newcodon)
       return tokens
       
      tokens.append(newcodon)
      pointer = pointer + CODONLENGTH
 
    return tokens
    
  def parse_tokens(self,tokenstring):
    """Finally a little semantics. Goes through tokens and assigns operators et al. to them.
       Numbers are put where numbers should be, fragments are Noneified, binding sites marked.
       Meaningless codons get "????" (from codon_to_operator)
       """
    parsedstring = copy.deepcopy(tokenstring)
    numberfollows = (self.operator_to_codon("META"),self.operator_to_codon("CATA"),self.operator_to_codon("SENS"),self.operator_to_codon("COPH"),self.operator_to_codon("JESU"))
    offcodons = (self.operator_to_codon("BIND"))
    
    pointer = 0
    while pointer < len(parsedstring):
      #ignore fragments
      if len(parsedstring[pointer]) < CODONLENGTH:
        parsedstring[pointer] = None
      #binding sites
      elif parsedstring[pointer] in offcodons:
        parsedstring[pointer] = self.codon_to_operator(parsedstring[pointer])
        pointer = pointer + 1
        parsedstring[pointer] = "*" + parsedstring[pointer].lower()     
      #numerical codes after specified codons
      elif parsedstring[pointer] in numberfollows:
        parsedstring[pointer] = self.codon_to_operator(parsedstring[pointer])
        pointer = pointer + 1
        parsedstring[pointer] = self.codon_to_integer(parsedstring[pointer])
      #normal codons
      else:
        parsedstring[pointer] = self.codon_to_operator(parsedstring[pointer])
      pointer = pointer + 1
    return parsedstring
  
  def text_to_gene(self,text):
    """converts a string into a gene representation. Useful for such as
     meaningful binding sequences"""
    message = ""
    for letter in text:
      message = message + self.integer_to_codon(ord(letter))      
    return message

  def gene_to_text(self,gene):
    """Reverses gene_to_text. Issues with non-letter characters."""
    plaintext = ""
    cdns = []
    pointer = 0
    while pointer < len(gene):
      cdn = ""
      subpointer = 0
      while subpointer < CODONLENGTH:
        cdn = cdn + gene[pointer+subpointer]
        subpointer = subpointer + 1 
      pointer = pointer + CODONLENGTH
      cdns.append(cdn)
    cdns = [chr(self.codon_to_integer(cd)+64) for cd in cdns]
    plaintext = "".join(cdns)
    return plaintext
            
class Genome(Code):
  def __init__(self,okay=[]):
    self.genome = okay
    self.tokenized = []
    self.parsed = []
    if okay:
      self.tokenize_genome()
      self.parse_genome()
      self.id_bindingsites()

  def id_bindingsites(self):
    bsites = []
    pointer = 0
    while pointer < len(self.parsed):
      if str(self.parsed[pointer])[0] == "*":
        bsites.append( (self.parsed[pointer],pointer) )
      pointer = pointer + 1    
    self.bindingsites = bsites
    return 1
    
  #=============================
  #GENOME CONSTRUCTING FUNCTIONS
  #=============================
     
  def add_gene(self,genestring):
    """Adds an arbitrary sequence to the end of the genome"""
    self.genome.extend( list(genestring.upper()) )
    return 1
    
  #probably obsolete because of add_operators
  def add_operator(self,operator):
    """Adds a specific codon to the end of the genome"""
    genestring = self.operator_to_codon(operator)
    self.add_gene(genestring)
    return 1

  def add_operators(self,codons):
    """Takes a tuple/list of operators and adds them in order"""
    for c in codons:
      self.add_operator(c)
    return 1
   
  
  #combine these two into one function eventually probably
  def add_random_genes(self,length):
    """Takes number of codons requested and adds them directly to genome.
       Generated randomly over whole alphabet."""
    randgen = []
    for i in xrange(length*CODONLENGTH):
      randgen.append(random.choice(basealphabet))
    randgen = "".join(randgen)
    self.add_gene(randgen)
    return 1

  def add_random_operators(self,length):
    """Takes number of codons requested and adds them directly to genome.
       Taken from operator list"""
    for i in xrange(length):
      thing = random.choice(operatorlist)
      genome. add_operator(thing)
    return 1
       
  #=================
  #PARSING FUNCTIONS
  #=================
  
  def tokenize_genome(self):
    self.tokenized =  self.tokenize_dna(self.genome)
    return 1
 
  def parse_genome(self):
    self.parsed = self.parse_tokens(self.tokenized)
    return 1
 
 
class Protein():
  def __init__(self,struct,phos=0,nom="",index=0):
    
    #static elements
    randnum = random.randrange(100)
    if nom:
      self.name = nom
    else:
      balls = random.choice(string.uppercase)+random.choice(string.lowercase)+random.choice(string.lowercase)
      appendletter = chr(65 + index)
      self.name = balls + appendletter
    
    self.id = randnum    
    self.structure = struct
    self.bindingsites = [operator for operator in self.structure if str(operator)[0] == "*"]

    #dynamic elements
    self.boundto = []
    self.phos = 0
  
  def phosphorylate(self, value=1):
    """alter phosphorylation status of protein by value"""
    self.phos = self.phos + value
    return self.phos
    
class Cell(Protein,Genome):
  def __init__(self,primer=["NULL"],memorysize=64,looplimit=4):
    
    #initialise the chemical array
    self.chemarray = []
    for register in range(memorysize):
      self.chemarray.append(0)

    #initialise the protein array
    self.protarray = []
    
    #add the Kadmon protein
    self.protkadmon = Protein(primer,0,"KadA")
    self.protarray.append(self.protkadmon)
    
    print "Priming with %s (Kadmon)! %s" % (self.protkadmon.name,primer)
    if primer[-1] != "LYSE":
      print "    CAUTION! %s (Kadmon) is not LYSE terminated!" % self.protkadmon.name
    print "Starting to execute %s (Kadmon)" % self.protkadmon.name
    self.execute_protein(self.protkadmon)
    print "Finished executing %s (Kadmon)" % self.protkadmon.name

    #if a protein without a LYSE is added, this stops it after looplimit terms.
    #Repress loop-limiting by setting looplimit = -1 upon instantiation
    
    lastprotein = self.protkadmon
    loopcaution = 0
    
    while self.protarray:
     currentprotein = self.protarray[0]
     #check to see if the same protein is at the top
     if currentprotein == lastprotein:
       #if so, increment the loop counter and warn
       loopcaution = loopcaution + 1
       print "%s executes again (%s times)" % (currentprotein.name, loopcaution)
     print "Starting to execute %s." % currentprotein.name
     self.execute_protein(currentprotein)
     print "Finished executing %s." % currentprotein.name
     lastprotein = currentprotein
     #halt unruly proteins
     if currentprotein == lastprotein:
       if loopcaution >= looplimit:
         print "Loop limit reached. Lysing %s externally." % currentprotein.name
         self.protarray.pop(0)
    
    multiprint("Tableau des proteines est vide. Fin.","Protein array empty. Stopping.")
     
    #display final chem memory state
    print ""
    print "Final chemical memory state:"
    if len(self.chemarray) != 64:
      print self.chemarray
    else:
      for i in xrange(8):
        if suppresstext:
          if max(self.chemarray) < 10:
            print "%s %s %s %s %s %s %s %s" % (self.chemarray[i*8],self.chemarray[i*8+1],self.chemarray[i*8+2],self.chemarray[i*8+3],self.chemarray[i*8+4],self.chemarray[i*8+5],self.chemarray[i*8+6],self.chemarray[i*8+7])
          elif max(self.chemarray) > 11 and max(self.chemarray) < 100:
            print "% 2s % 2s % 2s % 2s % 2s % 2s % 2s % 2s" % (self.chemarray[i*8],self.chemarray[i*8+1],self.chemarray[i*8+2],self.chemarray[i*8+3],self.chemarray[i*8+4],self.chemarray[i*8+5],self.chemarray[i*8+6],self.chemarray[i*8+7])
          else:
            print "% 3s % 3s % 3s % 3s % 3s % 3s % 3s % 3s" % (self.chemarray[i*8],self.chemarray[i*8+1],self.chemarray[i*8+2],self.chemarray[i*8+3],self.chemarray[i*8+4],self.chemarray[i*8+5],self.chemarray[i*8+6],self.chemarray[i*8+7])
        else:
          if max(self.chemarray) < 10:
            print "%s %s %s %s %s %s %s %s        %s %s %s %s %s %s %s %s" % (self.chemarray[i*8],self.chemarray[i*8+1],self.chemarray[i*8+2],self.chemarray[i*8+3],self.chemarray[i*8+4],self.chemarray[i*8+5],self.chemarray[i*8+6],self.chemarray[i*8+7],chr(self.chemarray[i*8]+64),chr(self.chemarray[i*8+1]+64),chr(self.chemarray[i*+2]+64),chr(self.chemarray[i*8+3]+64),chr(self.chemarray[i*8+4]+64),chr(self.chemarray[i*8+5]+64),chr(self.chemarray[i*8+6]+64),chr(self.chemarray[i*8+7]+64))
          elif max(self.chemarray) > 11 and max(self.chemarray) < 100:
            print "% 2s % 2s % 2s % 2s % 2s % 2s % 2s % 2s        %s %s %s %s %s %s %s %s" % (self.chemarray[i*8],self.chemarray[i*8+1],self.chemarray[i*8+2],self.chemarray[i*8+3],self.chemarray[i*8+4],self.chemarray[i*8+5],self.chemarray[i*8+6],self.chemarray[i*8+7],chr(self.chemarray[i*8]+64),chr(self.chemarray[i*8+1]+64),chr(self.chemarray[i*+2]+64),chr(self.chemarray[i*8+3]+64),chr(self.chemarray[i*8+4]+64),chr(self.chemarray[i*8+5]+64),chr(self.chemarray[i*8+6]+64),chr(self.chemarray[i*8+7]+64))
          else:
            print "% 3s % 3s % 3s % 3s % 3s % 3s % 3s % 3s        %s %s %s %s %s %s %s %s" % (self.chemarray[i*8],self.chemarray[i*8+1],self.chemarray[i*8+2],self.chemarray[i*8+3],self.chemarray[i*8+4],self.chemarray[i*8+5],self.chemarray[i*8+6],self.chemarray[i*8+7],chr(self.chemarray[i*8]+64),chr(self.chemarray[i*8+1]+64),chr(self.chemarray[i*+2]+64),chr(self.chemarray[i*8+3]+64),chr(self.chemarray[i*8+4]+64),chr(self.chemarray[i*8+5]+64),chr(self.chemarray[i*8+6]+64),chr(self.chemarray[i*8+7]+64))
  def change_chemarray(self,index,quantity):
    self.chemarray[index] = self.chemarray[index] + quantity
   
  def display_protarray(self):
    nice = ""
    for pro in self.protarray:
      nice = nice + " " + pro.name
    print "Protein Array: " + nice
    
  def execute_protein(self,prot):
    passcodons = ["NULL","FOOO","???","STOP","NULL","BIND","BEND","BEGN"]
    pointer = 0
    """Walks through a protein and executes it """
    while pointer < len(prot.structure):
      
      #PASS THROUGH PASSING CODONS IN passcodons      
      if prot.structure[pointer] in passcodons:
        multiprint("    %s fait rien (%s)" % (prot.name, prot.structure[pointer]),"    %s does nothing (%s)" % (prot.name, prot.structure[pointer]))
        pass

      #AUTOLYSIS!
      if prot.structure[pointer] == "LYSE":
        #lyse bound protein
        if prot.boundto:
          self.protarray.remove(prot.boundto) 
          print "    %s lyses bound %s." % (prot.name,prot.boundto.name)
          prot.boundto = []
        #or else itself
        else:
          self.protarray.remove(prot)
          multiprint("    %s se lyse." % prot.name,"    %s lyses itself." % prot.name)
          break

      #MAKING NEW PROTEINS
      if prot.structure[pointer] == "TRAN":
        
        #look for next binding site on protein
        nextbindingsite = ""
        pointer = pointer + 1
        while str(prot.structure[pointer])[0] != "*":
          pointer = pointer + 1
        nextbindingsite = prot.structure[pointer]
        
        #find matching binding site(s) on genome:
        relevantbindindices = []
        for entry in genome.bindingsites:
          if entry[0] == nextbindingsite:
            relevantbindindices.append(entry[1])
        
        multiprint("    %s trouve %s \"%s\" dans genome." % (prot.name, len(relevantbindindices), nextbindingsite),"    %s finds %s \"%s\" in genome." % (prot.name, len(relevantbindindices), nextbindingsite))
        
        for entry in relevantbindindices:
          transpointer = entry
          while genome.parsed[transpointer] != "BEGN":
            transpointer = transpointer + 1
          tempnewstruct = []
          while genome.parsed[transpointer] != "STOP":
            transpointer = transpointer + 1
            tempnewstruct.append(genome.parsed[transpointer])
          
          if not tempnewstruct:
            #fall back on this
            newstruct = ["NULL"]
          else:
            #remove trailing STOP as well
            newstruct = tempnewstruct[:-1]
          #make the new protein
          self.protarray.append(Protein(newstruct,0))
          #bind the new protein to the TRANsing protein
          prot.boundto = self.protarray[-1]

          multiprint("    Nouvelle Proteine %s! %s (offset du genome %s)" % (prot.boundto.name,prot.boundto.structure,entry),"    New Protein %s! %s (genome offset %s)" % (prot.boundto.name,prot.boundto.structure,entry))
          if prot.boundto.structure[-1] != "LYSE":
            multiprint("    ATTENTION! %s n'est pas n'est pas resilie par LYSE!" % prot.boundto.name,"    CAUTION! %s is not LYSE terminated!" % prot.boundto.name)
      #BACK STUFF NEEDS DOING PROPERLY
      if prot.structure[pointer] == "BACK":
       #unbind bound protein
       if prot.boundto:
         temp = prot.boundto.name
         prot.boundto = []
         multiprint("    %s detache %s!" % (prot.name,temp),"    %s unbinds %s!" % (prot.name,temp))
       else:
         print("    %s unbound nothing (BACK)" % (prot.name))
       #check to see if next codon is binding sequence
       if prot.structure[pointer+1] == "BIND":
         pointer = pointer + 2
         backbind = prot.structure[pointer]
         print "    %s activates binding site %s" % (prot.name,backbind)

         #remove self from the list of proteins under consideration
         goodlist = [protein for protein in self.protarray if protein.name != prot.name]
         #filter this list for proteins containing the binding site ID'd in backbind
         betterlist = [bro for bro in goodlist if backbind in bro.bindingsites]
         #bind to the first member of this array
         if betterlist:
           prot.boundto = betterlist[0]
           print "    %s binds %s through %s" % (prot.name,prot.boundto.name,backbind)
         #skip into BEND, the incrementer at the end will pick it up
         else:
           print "    Binding site %s on %s finds no target." % (backbind, prot.name)
         pointer = pointer + 1
       
      #PHOSPHORYLATE
      if prot.structure[pointer] == "PHOS":
        if prot.boundto:
          prot.boundto.phosphorylate()
          print "    %s increases phosphorylation of %s to %s" % (prot.name,prot.boundto.name,prot.boundto.phos)
        else:
          prot.phosphorylate()
          multiprint("    %s s'augment sa phosphorylation de %s" % (prot.name,prot.phos),"    %s increased phosphorylation to %s" % (prot.name,prot.phos))

      #DEPHOSPHORYLATE
      if prot.structure[pointer] == "DEPH":
        if prot.boundto:
          prot.boundto.phosphorylate(-1)
          print "    %s decreases phosphorylation of %s to %s" % (prot.name,prot.boundto.name,prot.boundto.phos)
        else:
          prot.phosphorylate(-1)
          multiprint("    %s se diminue sa phosphorylation de %s" % (prot.name,prot.phos),\
        "    %s decreased phosphorylation to %s" % (prot.name,prot.phos))

      #CONDITIONAL PHOSPHORYLATE
      if prot.structure[pointer] == "COPH":
        pointer = pointer + 1
        if self.chemarray[prot.structure[pointer]] != 0:
          if prot.boundto:
            prot.boundto.phosphorylate()
            print("    %s increased phosphorylation of %s to %s -- Register %s nonzero" % (prot.name,prot.boundto.name,prot.boundto.phos,prot.structure[pointer]))
          else:
            prot.phosphorylate()
            print("    %s increased phosphorylation to %s -- Register %s nonzero"\
           % (prot.name,prot.phos,prot.structure[pointer]) )
        else:
          print("    %s did nothing -- Register %s is zero" % (prot.name,prot.structure[pointer]))
          
      #METABOLISE A CHEMICAL
      if prot.structure[pointer] == "META":
        pointer = pointer + 1
        index = int(prot.structure[pointer])
        self.change_chemarray(index,prot.phos)
        multiprint("    %s augment le registre %s par %s." % (prot.name,index,prot.phos),\
        "    %s increased register %s by %s." % (prot.name,index,prot.phos))
      #CATABOLISE A CHEMICAL
      if prot.structure[pointer] == "CATA":
        pointer = pointer + 1
        index = int(prot.structure[pointer])
        self.change_chemarray(index,-prot.phos)
        multiprint("    %s diminue le registre %s par %s." % (prot.name,index,prot.phos),\
        "    %s decreased register %s by %s." % (prot.name,index,prot.phos))
      #PHOSPHORYLATE SELF ACCORDING TO CHEMICAL CONCENTRATION
      if prot.structure[pointer] == "SENS":
        pointer = pointer + 1
        if prot.boundto:
          prot.boundto.phos = self.chemarray[prot.structure[pointer]]
          print "    %s's phosphorylation is %s (from register %s by %s)" % (prot.boundto.name,prot.boundto.phos,prot.structure[pointer],prot.name)
        else:
          prot.phos = self.chemarray[prot.structure[pointer]]
          multiprint("    La phosphorylation de %s est %s (au registre %s)" % (prot.name,prot.phos,prot.structure[pointer]),\
        "    %s's phosphorylation is %s (from register %s)" % (prot.name,prot.phos,prot.structure[pointer]))        
      
      #SET CHEMICAL CONCENTRATION ACCORDING TO PHOSPHORYLATION STATE
      #"JESU n" : n -> register[phosphorylation]
      if prot.structure[pointer] == "JESU":
        pointer = pointer + 1
        if prot.boundto:
          self.chemarray[abs(prot.boundto.phos)] = prot.structure[pointer]
          print "    %s phosphor-sets register %s by %s to %s" % (prot.name,prot.boundto.phos,prot.structure[pointer],prot.boundto.name)
        else:
          self.chemarray[abs(prot.phos)] = prot.structure[pointer]
          print "    %s phosphor-sets register %s to %s" % (prot.name,prot.phos,prot.structure[pointer])
      #AND MOVE FORWARD
      pointer = pointer + 1
                
if __name__ == "__main__":
 
  #utility functions
  
  def prettyintro():
    """ so pretty :3 """
    print "       __  __  ___  _   _  ___  ____  "  
    print "      |  \/  |/ _ \| \ | |/ _ \|  _ \ "
    print "      | |\/| | | | |  \| | | | | | | |"
    print "      | |  | | |_| | |\  | |_| | |_| |"
    print "      |_|  |_|\___/|_| \_|\___/|____/ "
    print "      robert harry nicodemus williams~"
    print ""
    multiprint("                 Bonjour.","                 Hello.")
    print ""
    
  def debug():
    """Gets run when debug flag is added -- this just allows low-access level to the 
    contents of the genome and Kadmon primer at the level of code. """
    genome.add_operators(("NULL", "PHOS", "PHOS", "NULL", "BIND"))
    genome.add_gene("TATGAG")
    genome.add_operators(("BEND", "BEGN", "PHOS", "PHOS", "PHOS", "META", "FOOO", "FOOO","LYSE","STOP","NULL","BIND"))
    genome.add_gene("TATGAG")
    genome.add_operators(("BEND","BEGN","PHOS","META","PHOS","DEPH","META","FOOO","DEPH","NULL","BACK","BIND"))
    genome.add_gene("CATGAT")
    #genome.add_random_genes(2000000)
    #genome.add_operator("BEND")
    genome.add_operators(("BEND","PHOS","PHOS","PHOS","PHOS","DEPH","FOOO"))
    genome.add_operators(("BACK","LYSE","STOP"))
    genome.add_operators(("NULL","JESU","PHOS","BIND"))
    genome.add_gene("TATGAG")
    genome.add_operators(("BEND","BEGN","BIND"))
    genome.add_gene("CATGAT")
    genome.add_operators(("BEND","PHOS","PHOS","PHOS","PHOS","PHOS","PHOS","META"))
    genome.add_gene("AGA")
    genome.add_operators(("LYSE","STOP"))
    genome.add_operators(("NULL","BIND"))
    genome.add_gene("TATGAG")
    genome.add_operators(("BEND","BEGN","PHOS","PHOS","META"))
    genome.add_gene("CTA")
    genome.add_operators(("NULL","BIND"))
    genome.add_gene("ATTTTTTTT")
    genome.add_operators(("BEND","DEPH","LYSE","STOP"))
    
    genome.tokenize_genome()
    genome.parse_genome()
    genome.id_bindingsites()
    
    print "".join(genome.genome)
    cell = Cell(["TRAN", "BIND", "*tatgag", "BEND", "BACK", "NULL", "LYSE"])
    #cell = Cell(["PHOS", "PHOS", "PHOS", "META", 5 , "NULL", "SENS", 5, "META", 6, "LYSE"])
    
  def readdna(filename):
    """Reads a .dna file and returns a DNA string. Argument is the name of the file (with .dna)"""
    try:
      genomefile = open(filename,"r")
    except IOError:
      multiprint("%s n'existe pas!" % filename, "%s does not exist!" % filename)
      quit()  
    rawgenomecode = genomefile.readlines()
    genomefile.close()
    
    halfbaked = []
    for line in rawgenomecode:
      newline = [letter for letter in line if letter in "GACT"]
      halfbaked.extend(newline)
    genomecode = "".join(halfbaked)
    return genomecode
    
  def traduire():
    """Translates a .rna file into DNA string. Takes no arguments, instead asking for input and output filenames."""
    multiprint("Quel est le nom du fichier (avec ou sans .rna) qui contient les codons?","codon filename (with or without .rna)")
    codonfilename = raw_input("> ")      
    if codonfilename[-4:] != ".rna": codonfilename = codonfilename + ".rna"
    try:
      codonfile = open(codonfilename,"r")
    except IOError:
      multiprint("%s n'existe pas!" % codonfilename,"%s does not exist!" % codonfilename)
      quit()
    rawcodons = codonfile.read()
    codonfile.close()
    codons = rawcodons.split(" ")
    if codons[-1] == "": codons = codons[:-1]
    
    #convert codons back into dna
    dnastring = ""
    print codons
    for codon in codons:
      #normal case
      if codon in operatorlist:
        dnastring = dnastring+(code.operator_to_codon(codon))
      #binding sequence
      elif codon[0] == "*":
        foo = codon[1:].upper()
        dnastring = dnastring + foo
      #numerical codon
      elif codon[0] in string.digits:
        foo = code.integer_to_codon(codon)
        dnastring = dnastring + foo
      #revert to ballscodon
      else:
        dnastring = dnastring+random.choice(ballscodons)
    multiprint("Quelle est le nom (avec ou sans .dna) du fichier de sortie? Laisser vide pour le meme nom.","output filename? (with or without .rna) leave empty for same name")
    boules = raw_input("> ")
    writefilename = codonfilename[:-4]+".dna"
    if boules:
      if boules[-4:] != ".dna":
        writefilename = boules + ".dna"
    
    writefile = open(writefilename,"w")
    writefile.write(dnastring)
    writefile.close()
    quit()
  
  def numbersplease():
    """Prints out a list of numbers and their corresponding codons."""
    for num in xrange(len(basealphabet)**CODONLENGTH):
      print "%s %s" % (num,code.integer_to_codon(num))
    quit()
  
  def translatetext():
    """Allows user to input a string and recieve DNA string representing it."""
    print "Input plaintext."
    plaintext = raw_input("> ")
    print "DNA string for \"%s\":\n%s" % (plaintext,code.text_to_gene(plaintext))
    quit()
    
  def multiprint(frstring,enstring):
    """Multilingual (for multi in {french, english, caveman}) version of print(). Takes two arguments - the french and english
       versions of the same string. Caveman is emergent."""
    if language == "french":
      print frstring
    elif language == "english":
      print enstring
    else:
      #revert to caveman
      print "Ooga booga bekos."
  
  def complicatedinput():
    """old school type input mechanism, holds your hand a lot"""
    multiprint("Quel est le nom du fichier (avec ou sans .dna) qui contient le code source?","sourcecode filename? (with or without .dna")
    filename = raw_input("> ")      
    if filename[-4:] != ".dna": filename = filename + ".dna"
    
    #actually put contents of file into genome
    genome = Genome(readdna(filename))
    
    #write .rna file
    rnafile = open("%s.rna"%filename[:-4],"w")
    for codon in genome.parsed:
      rnafile.write("%s " % str(codon))
    rnafile.close()
    
    #and now to compute we need the kadmon primer
    multiprint("Quel est le nom du fichier (avec ou sans .rna) qui contient l'amorce?","primer filename? (with or without .rna")
    primerfilename = raw_input("> ")
    if primerfilename[-4:] != ".rna": primerfilename = primerfilename + ".rna"
    primerfile = open(primerfilename,"r")
    rawprimer = primerfile.read()
    primerfile.close()
    primer = rawprimer.split(" ")[:-1]

    # all the input we need has been given and we can finally set the cell in its inexorable motion
    cell = Cell(primer)
    sys.exit(0)
  
  def displayhelp():
    """displays etc..."""
    multiprint("nom de dieu jsais pas cest embrouille ne me demandez pas jne travaille quici","o god i dont know its complicated dont ask me i just work here")
    quit()
  #======================================================================
  #START MAIN PROPER
  #======================================================================
  
  global language
  language = "french"
  
  #whether to print the text in the memorydump or no, 1 = yes, 0 = no
  global suppresstext
  suppresstext = 0
  
  if "-s" in sys.argv:
    suppresstext = 1
    
  #initialise the code, which is always useful
  code = Code()

  #USE ENGLISH YOU FUCKING IMMIGRANTS
  if "-e" in sys.argv:
    language = "english"

  #display the nice title in french, english or caveman
  prettyintro()
  
  #O GOD HALP HALP
  if "-h" in sys.argv:
    displayhelp()
      
  #.rna => .dna mode
  if "-t" in sys.argv:
    traduire()
  
  #print numberlist
  if "-n" in sys.argv:
    numbersplease()
  
  #quick version of -f
  if "-x" in sys.argv:
    xat = sys.argv.index("-x")
    arguments = [sys.argv[xat+1],sys.argv[xat+2]]
    files =  arguments[-2:]
    rnafilename = [fil for fil in files if fil[-3] == "r"][0]
    dnafilename = [fil for fil in files if fil[-3] == "d"][0]

    genome = Genome(readdna(dnafilename))
    
    primerfile = open(rnafilename,"r")
    rawprimer = primerfile.read()
    primerfile.close()
    primer = rawprimer.split(" ")[:-1]
    
    cell = Cell(primer)
    
    sys.exit(0)
  
  #translate text into dna string
  if "-c" in sys.argv:
    translatetext()
    
  #.dna => .rna => computation!
  if "-f" in sys.argv:
    complicatedinput()
    
  # debug by default - ie. what i've been using up to now...
  else:
      genome = Genome()
      debug()
