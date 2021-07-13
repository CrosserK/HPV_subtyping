# Info
FindHPVMutations er et script der kan finde mutationer i HPV typer på nukleotid og aminosyre niveau. Der inputtes hvilke HPV typer der er i en prøve (f.eks gennem info fra en coverage analyse) og scriptet vil så finde varianter for hver type. Scriptet kan også finde hvilken HPV type der er tilstede, men det anbefales ikke hvis der er flere HPV typer i en prøve. 

# Fremgangsmåde
1. Læg alle fastq filer fra ION Torrent ind i mappen /home/pato/Skrivebord/HPV16_projekt/FASTQ
Der kan godt oprettes undermapper i mappen for at holde styr på andre FASTQ filer - scriptet søger ikke i undermapper. 
Hvis der har været tidligere ufuldstændige kørsler eller enkelte moduler har været kørt, skal alle eventuelle filer der ender med _filt.fastq fjernes fra denne mappe. Åben terminal (ctrl+alt+t, og indsæt følgende for at slette dem: 
rm -I /home/pato/Skrivebord/HPV16_projekt/FASTQ/*_filt.fastq; 
og tryk enter. Den spørger om du vil slette x antal filer, skriv y og tryk enter hvis ja (y=yes)). 

2. Åben tekstfilen /home/pato/Skrivebord/HPV16_projekt/FindHPVMutations_options.txt.
I denne tekstfil kan indstillinger for kørslen sættes. Vær opmærksom på ikke at ændre i navnene til venstre for ”=” tegn. Se Indstillinger sektion for forklaringer af de forskellige muligheder. 

3. Som udgangspunkt angives hvilke HPVtyper der er i en prøve, angives disse for hver fastqfil i en individuel tekstfil i /References/InputRefs/<fastqnavn>.txt. Hver type skal stå på en linje for sig. Undgå mellemrum på endelser. Dette er dog ikke nødvendigt hvis der køres hvis der køres VirStrain subtypering.
4. Når ovenstående parametre er sat , gemmes den ændrede FindHPVMutations_options.txt fil (ctrl+s). Tekstfilen FindHPVMutations_options.txt kan nu lukkes. 
  
5. Åben terminalen (ctrl+alt+t) og start programmet ved at skrive: FindHPVMutations.sh og trykke enter. 
  

# Indstillinger
**RunName** er et navn man kan give hele kørslen. Kørslen kommer til at hedde <RunName>_<klokkeslæt>_<dato>
  
**MainF** er hovedmappen, hvorfa alt køres og alle undermapper ligger. Scriptet er ikke ligenu flytbart, men få ting skal fikses før det kan blive det. 
  
**Qualtrim** er en kvalitetsscore som cutadapt programmet vil bruge til at fjerne baser fra 3' enden af reads, hvis de har en score under den angivne værdi
  
**MinLen** er min. længde reads skal være, ellers bliver de filtreret fra
  
**MaxLen** er max. længde reads skal være, ellers bliver de filtreret fra
  
**AmpliconRef** er den reference amplicons er lavet ud fra og bruges til at kunne klippe regioner af reads udenfor amplicon væk. 
  
BedFileNameX er navn på hypotetisk bedfil med komplimentære regioner til at cutte primer sekvenser væk fra amplicons. Der skal ligges en fil med endelsen compl1 og en med endelsen compl2 der udgør bed fil for hver amplicon primer pool. Det angivne navn skal da ende med _compl, så finder scriptet det to bed filer med endelserne 1 og 2. 
  

De to sidstenævnte variabler bliver ikke brugt medmindre der bliver alignet til den eksakte angivne ampliconreference. 


## Modul indstillinger
**customRefForAll** [true/false] Hvis true ignoreres inputrefs og VirStrain subtypering og alle fastsqfiler bliver alignet til samme reference(r) angivet i cRef. Slå modulerne CombineRefs, SplitAndVarCall og VirStrainGenoAndSubTyping fra, hvis denne aktiveres.
  
**cRef** Her angives hvilke reference(r) der skal alignes til, hvis customRefForAll er aktiveret. Flere referencer angives med et tab i mellem, eks: 
HPV35_X74477_1	HPV70_U21941_1	HPV49_X74480_1
  
**customName** [true/false] Hvis true, ignoreres RunName og et angivet navn i cName bruges i stedet, hvor der ikke sættes dato på. Fordelen er at det kan bruges til at køre specifikke moduler på en evt. ældre kørsel, der allerede har et navn med klokkeslet
  
**cName** se ovenstående.
  
**indexReferences** [true/false] hvis true indekseres alle referencer som ligger i /References hvis der ikke allerede findes en .fasta fil for referencen i References/IndexedRef. Fasta kopieres til undermappe af samme navn i References/IndexedRef og indeks filer genereres i samme mappe.
  
**cutOutsideAmplicons** [true/false] Klipper alle sekvenser uden for amplicon bedfil. Man vil risiskere at miste data hvis prøve har HPV type andre end den brugt til at lave ampliconpanel
  
**qualityFilt** [true/false]  angiver om fastqfiler skal filtreres. Dette modul skal være aktiveret, hvis der skal alignes, da align og variantcall script kun leder efter filtrerede filer. 
  
**CombineRefs** [true/false] skal være true hvis der er angivet flere genotyper i en eller flere fastqfiler
  
**AlignAndVarCall** [true/false] modul der aligner og laver variantcalls
  
**SplitAndVarCall** [true/false] modul der splitter og laver nye variantcalls på individuelle fastqfiler der indeholder flere angivne HPV typer 
  
**RunAnnoRSCript** [true/false] Modul der annoterer variantcalls med aminosyre ændringer
  
**RunSiteCovScript**  [true/false] Modul der tæller dybden af reads på hver variant position 
  
**RunNoCallsScript**  [true/false]  Modul der sammenligner hvor mange af de brugte fastqfiler som er alignet til en given reference der ikke har tilstrækkelig dybde på variant position. Dette modul køres til allersidst og giver desuden en summerende tabel med resultater for alle varianter og referencer. 
  

**VirStrainGenoAndSubTyping** [true/false] Angiver om VirStrain skal bruges til at finde hoved og subtyper før alignment og variantcalling i stedet for manuelt angivne typer. Ikke optimalt hvis der er flere typer i en prøve. Modulerne CombineRefs og SplitAndVarCall skal deaktiveres hvis denne er aktiveret.
  
**VirStrainMaindb** Navn på VirStrain genereret database over hovedtyper. Skal ligge i /References/
  
**VirStrainSubFolders** Navn på mappe der indeholder VirStrain genereret subtypedatabaser i individuelle mapper. Skal ligge i /References/ og database for hver undertype skal ligge i en mappe for sig under denne mappe. 
  
## Avancerede indstillinger
### Alignment indstillinger
For at skrue på aligneren BWA’s parametre åbnes scriptet AlignAndVariantCall.sh i en tekst editor som wordpad eller sublime text (anbefales). Herefter scrolles ned til linje 85, som starter med teksten ”bwa mem”. Nye parametre tilføjes lige efter teksten bwa mem i linje 85, adskilt af mellemrum. Parameter indstillinger for bwa mem kan ses på http://bio-bwa.sourceforge.net/bwa.shtml under Commands and options → mem → Options.  
### Variant calling indstillinger
For at skrue på variant calleren GATK HaplotypeCallers parametre åbnes scriptet Scripts/Bash/AlignAndVariantCall.sh i en tekst editor som wordpad eller sublime text (anbefales). Herefter scrolles ned til linje 135 og 181, som starter med teksten ”gatk --java-options "-Xmx4g" HaplotypeCaller”. Her skal det noteres at linjen slutter med en backslash ”\” hvilket betyder at næste linje læses som en del af foregående. Nye parametre tilføjes til sidste linje af de \  sammenhængende linjer eller med en ny linje ved at sætte \ på sidste sammenhængende linje for at forbinde den til den nye linje. Det er vigtigt at man er opmærksom på at ændre de samme ting for både linje 135 og 181. Gem fil inden kørsel. Parameter indstillinger for variantcalleren kan ses på https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller
### Variant filtrering
For at skrue på variant calleren GATK HaplotypeCallers parametre åbnes scriptet Scripts/Bash/AlignAndVariantCalll.sh i en tekst editor som wordpad eller sublime text (anbefales). Herefter scrolles ned til linje 226 som starter med ”gatk --java-options "-Xmx4g" VariantFiltration”. Her kan ændres i parametre i linjerne under, som er forbundet til 286 med backslash ” \” ved at ændre i tallet efter < eller > tegnet. Gem fil inden kørsel. Hvert filter er forklaret på https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
### Nye subtyper med VirStrain
Hvis der skal subtyperes nye typer, skal en ny multiple sequence alignement efterfulgt af en generering af en database til VirStrain foretages. Se MAFFT&VirStrain.sh script i Scripts/Bash mappen. Herefter henvises VirStrain til de nylavede databaser ved at ændre i linjerne 
  
”VirStrainMaindb=”
  
”VirStrainSubFolders=”
  
i filen ”FindHPVMutations_options.txt”
  
Hvor den første er en database med alle HPV hovedtyper og den anden er en overmappe med undermapper af databaser med alle HPV undertyper. Hvis en ny undertype tilføjes til eks HPV16 og er lagt ind i den samme database under subtyper, som brugt ved tidligere kørsel, behøves de to ikke ændres. 
Individuelle referencer skal derefter lægges ind i /home/pato/Skrivebord/HPV16_projekt/References 
Det er nødvendigt at der ligeledes er en tilhørende gff3 filer til hver ny reference i /References/GFFfiles hvis der skal kaldes aminosyre ændringer. 
Sæt indexReferences=true i FindHPVMutations_options.txt, første gang der køres med nye referencer. Der bliver genereret nye index filer for hver reference og de bliver gemt i undermapper til fremtidige kørsler. 
  
# Resultater
Resultater kommer ud i filen Annotation_results/AnnotationFrequency_FASTQfiles_<RunName>_<time>_<date>_AllResults.txt, som indeholder summerende statistik over alle fundne varianter.
  
Filen Annotation_results/Annotations_FASTQfiles_<RunName>_<time>_<date>.txt indeholder info over hvilke varianter der er kaldt i hver enkelt fastqfil. 
  
Alignede reads ligger i bamfiler og varianter i  vcf filer som er lagt i  Results/<RunName>/<fastqnavn>/<reference>/
  
En summerende fil for antal varianter fundet i en fastqfil ligger i Results/<RunName>/<fastqnavn>/MismatchCounts_filt_<fastqnan>_run, hvor første kolonne er referencenavn og anden kolonne er antal fundne mismatch. Et lavt antal mismatch betyder ikke nødvendigvis at reads matcher reference perfekt, det kan også betyde at der ikke er alignet noget og der derfor ikke findes nogen mismatches.
  
Bam filer med endelsen ".sort.bam" har ikke fået markeret duplikat reads
  
Bam filer med endelsen ".sort.dup.readGroupFix.bam" har fået markeret duplikat reads
  
VCF filer med endelsen ".sort.dup.readGroupFix_filtered.vcf" har fået filtreret variant calls, men indeholder dem stadig, så man kan se hvorfor de er blevet filtreret.
  
VCF filer med endelsen ".sort.dup.readGroupFix_filtered_FiltEx_headerfix.vcf" har fået fjernet filtrerede variant calls helt. 
  
Hvis VirStrain subtypering er aktiveret ligger resultater for subtypering i VirStrain_run/<RunName>/VirStrain_summary.txt og i undermapper med samme navn som den enkelte fastqfil i VirStrain_report.txt. 
  
Resultater for genotypering for hovedmapper ligger på tilsvarende vis under mappen VirStrain_run/<RunName>/Genotypecalls.

Efter en kørsel er der lavet filer med navn FASTQfiles_<RunName>_<time>_<date>.txt og  FASTQfiles_<RunName>_<time>_<date>_runnames.txt i hovedmappen. Disse filer har været brugt til at finde fastqfiler som bruges i kørslen og efterlades i tilfælde af at der skal laves modificeringer eller enkelte moduler skal genkøres. De kan slettes manuelt efter en overstået kørsel.
