# encoding: utf-8
#!/usr/bin/env ruby

gem 'rest-client', '~> 2.0'
require 'rest-client'
require 'zlib'
require 'open-uri'
require 'json'

#procedura per inserire i dati nella hash nel caso di gene db
#riceve in input la hash, il codice univoco del trascritto, il cromosoma, numero esoni nel gene, numero trascritti, numero esoni nel gene, locazione, array posizioni degli esonie e se è una cds o meno 
#se il codice non ha .1, .2 o .3 aggiungo .1 per avere coerenza tra i codici dei due db
def inserimento (h, cod, c, neg, nt, net, p, arr, cdsb)
	if ((cod.include?(".1"))or(cod.include?(".2"))or(cod.include?(".3")))!=true then
		cod=cod+".1"
	end
	array=[c, neg, nt, net, p, arr, cdsb]
	h[cod]=array
end

def verificaesoni (primo, secondo)
	esonicomuni=0
	primo.each do |v|
		if secondo.include?(v) then
			esonicomuni=esonicomuni+1
		end					
	end
	esoni=esonicomuni+(primo.length-esonicomuni)+(secondo.length-esonicomuni)

	return esoni
end

#procedura per verificare il numero di esoni totale tra tre diverse isoforme
#riceve in input gli array con le posizioni dei diversi esoni nella sequenza
def verificaesonitre (primo, secondo, terzo)
	esonicomuni=0
	#se il secondo e il terzo includono lo stesso esone lo conto come comune
	primo.each do |v|
		if secondo.include?(v) then
			if terzo.include?(v) then 
				esonicomuni=esonicomuni+1
			end
		end					
	end
	primosecondo=0
	primo.each do |v|
		if secondo.include?(v) then
			if terzo.include?(v)==false then 
				primosecondo=primosecondo+1
			end
		end					
	end
	secondoterzo=0
	secondo.each do |v|
		if terzo.include?(v) then
			if primo.include?(v)==false then 
				secondoterzo=secondoterzo+1
			end
		end					
	end
	primoterzo=0
	primo.each do |v|
		if terzo.include?(v) then
			if secondo.include?(v)==false then 
				primoterzo=primoterzo+1
			end
		end					
	end
	#il numero di esoni totali è così calcolato
	esoni=esonicomuni+primosecondo+secondoterzo+primoterzo+(primo.length-esonicomuni-primosecondo-primoterzo)+(secondo.length-esonicomuni-primosecondo-secondoterzo)+(terzo.length-secondoterzo-esonicomuni-primoterzo)
	return esoni
end

#variabili per estrarre i nomi dei file nei db
codiciplasmodb=""
codicigenedb=""
primagenedb=/\suniqueName="/ #pattern prima del codice del file del cromosoma nell'elenco dei nomi dei file
dopogenedb=/"/
primaplasmo=/<record\sid='/ #pattern prima del codice del file del cromosoma nell'elenco dei nomi dei file
dopoplasmo=/'>/

#creo due array per contenere i dati dei 14 cromosomi
arrayplasmo= Array.new(14)
arraygene=Array.new(14)

#variabili per il pm in plasmodb
idarray=/source_id/
idarraybool=false
seqid=/sequence_id/
idseqbool=false
esoninelgene=/gene_exon_count/
esonigenebool=false
trascrittinelgene=/gene_transcript_count/
genetrasbool=false
countexon=/exon_count/
esonicont=false
location=/location_text/
posizione=false
cdata=/<!\[CDATA\[/
finecdata=/\]\]>/
sequence_id=""
esoniingene=""
trascrittigene=""
esoniintrasc=""
indice=""
locationgenomic=""

#variabili usate per raccogliere dati genedb
darivedere=""
contenutogz=""
codice=""
cdsjoin=/^FT\s+CDS\s+join\(/
parentesi=/\)/
locus=/^FT\s+\/locus_tag="/
other=/^FT\s+\/other_transcript="/
cds=/^FT\s+CDS\s+(\d+)\.\.(\d+)/
trna=/^FT\s+tRNA\s+(\d+)\.\.(\d+)/
trnacom=/^FT\s+tRNA\s+complement\((\d+)\.\.(\d+)/
cdscom=/^FT\s+CDS\s+complement\(/
cdscomplemjoin=/^FT\s+CDS\s+complement\(join\(/
doppiaparentesi=/\)\)/
ncrnacompl=/^FT\s+ncRNA\s+complement\((\d+)\.\.(\d+)/
ncrna=/^FT\s+ncRNA\s+(\d+)\.\.(\d+)/
rrna=/^FT\s+rRNA\s+(\d+)\.\.(\d+)/
rrnacomp=/^FT\s+rRNA\s+complement\((\d+)\.\.(\d+)/
cdsmot=/^FT\s+CDS_motif/
misc=/^FT\s+misc_feature/
rreg=/^FT\s+repeat_region/

#variabil usate nella ricerca dei casi particolari
intpos=/:/
strn=/\([\+-]\)/
casiparticolari=""
posinizioplasmo=""
posfineplasmo=""
posiniziogene=""
posfinegene=""
trasprob=""

#scarico da plasmodb i nomi dei file riguardanti p.falciparum
id = RestClient.get 'http://plasmodb.org/webservices/GenomicSequenceQuestions/SequencesByTaxon.xml?organism=Plasmodium%20falciparum%203D7'

#crea una stringa contenente tutti i nomi dei file
id.each_line do |linea|
	if linea=~primaplasmo then
		restante="#{$'}" #prendo la parte dopo
		if restante=~dopoplasmo then
			codice="#{$`}" #prendo la parte prima
			codiciplasmodb=codiciplasmodb+codice+"  "
		end
	end
end
			
#splitto la stringa creata, ottengo array contenente i codici
arrcodplas=codiciplasmodb.split("  ")

#scarico da genedb i nomi dei file riguardanti p.falciparum
uniqname= RestClient.get 'http://www.genedb.org/services/regions/inorganism.xml?organism=Pfalciparum'

#crea una stringa contenente tutti i nomi dei file
until uniqname==""
	if uniqname=~primagenedb then		
		restante="#{$'}"
		if restante=~dopogenedb then
			codice="#{$`}"
			uniqname="#{$'}"
			codicigenedb=codicigenedb+codice+"  "
		end
	else
		uniqname=""
	end
end
			
#splitto la stringa creata, ottengo array contenente i codici
arrcodgene=codicigenedb.split("  ")
				
#controllo corrispondenze
#corrispondenzacodici=true
#arrcodplas.each {|cod| if arrcodgene.include?(cod)==false then corrispondenzacodici=false end}
#if corrispondenzacodici==true
	#puts "i due db contengono gli stessi file"
#else 
	#puts "i due db non contengono gli stessi file"
#end

iniziorichiestamin="http://plasmodb.org/plasmo/webservices/GeneQuestions/GenesByLocation.xml?organismSinglePick=Plasmodium%20falciparum%203D7&chromosomeOptional=0"
iniziorichiesta="http://plasmodb.org/plasmo/webservices/GeneQuestions/GenesByLocation.xml?organismSinglePick=Plasmodium%20falciparum%203D7&chromosomeOptional="
finerichiesta="&start_point=1&end_point=0&o-fields=source_id,sequence_id,gene_exon_count,gene_transcript_count,exon_count,location_text"

#inserisco i dati dei 14 cromosomi sulla base dei file di plasmodb
1.upto(14) do |n|
	#creo hash nuova per ogni cromosoma
	hashplasmo= Hash.new
	if n<10 then
		nuovarichiesta= RestClient.get iniziorichiestamin+"#{n}"+finerichiesta
	else
		nuovarichiesta= RestClient.get iniziorichiesta+"#{n}"+finerichiesta
	end
	inserire=false
	#per ogni cromosoma
	nuovarichiesta.each_line do |linea|
		#se una variabile booleana è a true significa che nella riga successiva attendo quel contenuto
		if linea=~idarray then
			idarraybool=true
			inserire=false
		elsif linea=~seqid then
			idseqbool=true
		elsif linea=~esoninelgene then
			esonigenebool=true
		elsif linea=~trascrittinelgene then
			genetrasbool=true
		elsif linea=~countexon then
			esonicont=true
		elsif linea=~location then
			posizione=true		
		end
		#se la linea contiene dati e abbiamo una determinata var a true, cambiamo a false ed estraiamo il valore atteso
		#indice è il codice univoco del trascritto
		if(linea=~cdata) and (idarraybool==true) then
			idarraybool=false
			restante="#{$'}"
			if restante=~finecdata then
				indice="#{$`}"
			end
		#sequence_id è il cromosoma di appartenenza
		elsif (linea=~cdata) and (idseqbool==true) then
			idseqbool=false
			restante="#{$'}"
			if restante=~finecdata then
				sequence_id="#{$`}"										
			end
		elsif (linea=~cdata) and (esonigenebool==true) then
			esonigenebool=false
			restante="#{$'}"
			if restante=~finecdata then
				esoniingene="#{$`}"										
			end
		elsif (linea=~cdata) and (genetrasbool==true) then
			genetrasbool=false
			restante="#{$'}"
			if restante=~finecdata then
				trascrittigene="#{$`}"										
			end
		elsif (linea=~cdata) and (esonicont==true) then
			esonicont=false
			restante="#{$'}"
			if restante=~finecdata then
				esoniintrasc="#{$`}"										
			end
		elsif (linea=~cdata) and (posizione==true) then
			posizione=false
			restante="#{$'}"
			if restante=~finecdata then
				locationgenomic="#{$`}"										
			end
			#sappiamo he dopo aver travato la locazione il record è completo e può essere inserito
			inserire=true
		end							
		if inserire==true then
			#struttura del value [sequence_id, esoniingene, trascrittinelgene, esoni in trascritti, location]
			arraydati=[sequence_id, esoniingene.to_i, trascrittigene.to_i, esoniintrasc.to_i, locationgenomic]
			hashplasmo[indice]=arraydati						
		end
	end
	#inserisco nell'array di plasmodb la hash riguardante i dati del cromosoma
	arrayplasmo[n-1]=hashplasmo	
end
			
#download dati genedb
url="ftp://ftp.sanger.ac.uk/pub/project/pathogens/malaria2/3D7/3D7.latest_version/version3.1/2017/February_2017/"
estensione=".embl.gz"
c=0
arrcodgene.each do |cod| 
	download = open(url+cod+estensione)
	IO.copy_stream(download, cod+estensione)
	Zlib::GzipReader.open(cod+estensione) do |gz|
		loc=""
		str=""
		cdsbool=false
		a=[]
		aa=[]
		attiva=false #inizialmente a false per non inserire nulla prima di un trascritto valido 
		inserire=false
		motif=false #sarà true quando mi trovo in una cds_motif e non mi permetterà l'inserimento
		numeroesonigene=0
		numerotrascritti=0
		numeroesonitrascr=0
		hashgene= Hash.new
		contenutogz=gz.read
		#per ogni riga dell'embl verifico tramite pm con che situazione si ha a che fare per determinare e calcolare i vari valori e parametri del trascritto
		contenutogz.each_line do |linea|
			linea.chomp!
			if linea=~cdsjoin then
				attiva=true
				motif=false
				if inserire==true then
					inserimento(hashgene, codice, cod, numeroesonigene, numerotrascritti, numeroesonitrascr, loc, aa, cdsbool)
					numeroesonigene=0
					numerotrascritti=0
					numeroesonitrascr=0
					loc=""
					aa=[]
					a=[]
					inserire=false
					cdsbool=false
				end
				str="(+)"
				cdsbool=true
				restante="#{$'}"
				if restante=~parentesi then
					datiesoni="#{$`}"
					a=datiesoni.split(",")
					numeroesonitrascr=a.length
					a.each {|elemento| aa.push(elemento.split(".."))}
					loc=cod+":"+aa.first.first+".."+aa.last.last+str
				end
			elsif linea=~cdscomplemjoin then
				restante="#{$'}"
				attiva=true
				motif=false
				if inserire==true then
					inserimento(hashgene, codice, cod, numeroesonigene, numerotrascritti, numeroesonitrascr, loc,aa,cdsbool)
					numeroesonigene=0
					numerotrascritti=0
					numeroesonitrascr=0
					loc=""
					aa=[]
					a=[]
					inserire=false
					cdsbool=false
				end
				str="(-)"
				cdsbool=true
				if restante=~doppiaparentesi then
					datiesoni="#{$`}"
					a=datiesoni.split(",")
					numeroesonitrascr=a.length	
					a.each {|elemento| aa.push(elemento.split(".."))}
					loc=cod+":"+aa.first.first+".."+aa.last.last+str
				end
			elsif linea=~cdscom then
				motif=false
				attiva=true
				if inserire==true then
					inserimento(hashgene, codice, cod, numeroesonigene, numerotrascritti, numeroesonitrascr, loc, aa, cdsbool)
					numeroesonigene=0
					numerotrascritti=0
					numeroesonitrascr=0
					loc=""
					aa=[]
					a=[]
					cdsbool=false
					inserire=false
				end
				str="(-)"
				cdsbool=true
				restante="#{$'}"
				if restante=~parentesi then
					datiesoni="#{$`}"
					a=datiesoni.split(",")
					numeroesonitrascr=a.length
					aa=[]
					a.each {|elemento| aa.push(elemento.split(".."))}
					loc=cod+":"+aa.first.first+".."+aa.last.last+str
				end
			elsif (linea=~cds) then
				motif=false
				attiva=true
				if inserire==true then
					inserimento(hashgene, codice, cod, numeroesonigene, numerotrascritti, numeroesonitrascr, loc, aa, cdsbool)
					numeroesonigene=0
					numerotrascritti=0
					numeroesonitrascr=0
					loc=""
					aa=[]
					a=[]
					inserire=false
					cdsbool=false
				end
				str="(+)"
				cdsbool=true
				numeroesonigene=1
				numerotrascritti=1
				numeroesonitrascr=1
				loc=cod+":#{$1}..#{$2}"+str
				aa=["#{$1}","#{$2}"]
			elsif (linea=~ncrna) or (linea=~rrna) or (linea=~trna) then
				motif=false
				attiva=true
				if inserire==true then
					inserimento(hashgene, codice, cod, numeroesonigene, numerotrascritti, numeroesonitrascr, loc, aa,cdsbool)
					numeroesonigene=0
					numerotrascritti=0
					numeroesonitrascr=0
					loc=""
					aa=[]
					a=[]
					inserire=false
					cdsbool=false
				end
				str="(+)"
				numeroesonigene=1
				numerotrascritti=1
				numeroesonitrascr=1
				loc=cod+":.."+str
				aa=["#{$1}","#{$2}"]
			elsif (linea=~trnacom) or (linea=~ncrnacompl) or (linea=~rrnacomp) then
				motif=false
				attiva=true
				if inserire==true then
					inserimento(hashgene, codice, cod, numeroesonigene, numerotrascritti, numeroesonitrascr, loc, aa,cdsbool)
					numeroesonigene=0
					numerotrascritti=0
					numeroesonitrascr=0
					loc=""
					aa=[]
					a=[]
					inserire=false
					cdsbool=false
				end
				str="(-)"
				numeroesonigene=1
				numerotrascritti=1
				numeroesonitrascr=1
				loc=cod+":.."+str
				aa=["#{$1}","#{$2}"]
			elsif (linea=~cdsmot) or (linea=~misc) or (linea=~rreg) then
				motif=true
				if inserire==true then
					inserimento(hashgene, codice, cod, numeroesonigene, numerotrascritti, numeroesonitrascr, loc, aa,cdsbool)
					numeroesonigene=0
					numerotrascritti=0
					numeroesonitrascr=0
					loc=""
					aa=[]
					a=[]
					inserire=false
					cdsbool=false
				end
			elsif linea=~other then
				restante="#{$'}"
				if restante=~dopogenedb then
					altrotrascritto="#{$`}"
				end
				darivedere=darivedere+altrotrascritto+" "
			#se trovo un locus posso inserire i dati solo se non ho a che fare con una cds_motif
			elsif (linea=~locus) and (motif==false) then
				#locus contiene il nome univoco del trascritto
				#se attiva è true so di avere i dati necessari per poter inserire prima di sovrascrivere le vriabili
				if attiva==true then
					inserire=true
				end
				restante="#{$'}"
				if restante=~dopogenedb then
					codice="#{$`}"
				end
			end
		end
		#inserisco la hash del cromosoma
		arraygene[c]=hashgene
		#aggiorno il numero del cromosoma in considerazione
		c=c+1
	end
end

#gestione other_transcript
danonrifare=[]
chiavedue=""
chiaveuno=""
#per ogni trascritto versifico quante isoforme ha e in base a questo calcolo numero totale di esoni nel gene e trascritti nel gene
arraygene.each do |elem|
	elem.each do |key, value|
		if key.include?(".3") then
			#se ho tre isoforme le onsidero tutte e tre subito e non dovro riconsiderarle come coppia di isoforme
			chiavedue=key.sub(/\.3/, ".2")
			danonrifare.push(chiavedue)
			chiaveuno=key.sub(/\.3/, ".1")
			#numero di isoforme =3
			elem[chiavedue][2]=3
			elem[chiaveuno][2]=3
			elem[key][2]=3
			#verifico tramite una procedura esterna il numero di esoni nel gene
			esonig=verificaesonitre(elem[chiaveuno][5], elem[chiavedue][5], elem[key][5])
			elem[chiavedue][1]=esonig
			elem[chiaveuno][1]=esonig
			elem[key][1]=esonig
		end
		if key.include?(".2") then
			if danonrifare.include?(key)!=true then
				chiaveuno=key.sub(/\.2/, ".1")
				elem[chiaveuno][2]=2
				elem[key][2]=2
				esonig=verificaesoni(elem[chiaveuno][5], elem[key][5])
				elem[chiaveuno][1]=esonig
				elem[key][1]=esonig
			end
		end
	end
	#svuoto i trascritti da non rifare quando passo da un cromosoma all'altro
	danonrifare.clear
end
		
#per tutti i trascritti senza other_transcript stabilisco che il numero totale di esoni nel gene è uguali a quelli del trascritto in quanto è unica isoforma		
arraygene.each do |elem|
	elem.each do |key, value|
		if (value[1]==0) and (value[2]==0) then
			value[1]=value[3]
			value[2]=1
		end
	end
end

#stampiam su file tutti i trascriti derivanti da cds
stampa=""
arraygene.each do |elem|
	elem.each do |key, value|
		if (value[6]==true)
			stampa=stampa+key+" \n"
		end
	end
end
File.open ARGV[0], "w" do |output|
    output << stampa	
end

#cerchiamo casi particolari di non coincidenza tra le informazioni date dai due db
0.upto(13) do |pox|
	p=arrayplasmo[pox]
	g=arraygene[pox]
	p.each do |k,v|
		if g.has_key?(k) then
			sottoarr=g[k][4]
			value=v[4]
			trasprob=g[k][2]
			if sottoarr!=value then
				if value=~intpos then		
					restante="#{$'}"
					posizioniarr=restante.split("..")
					posinizioplasmo=posizioniarr[0]
					if posizioniarr[1]=~strn then	
						posfineplasmo="#{$`}"
					end
				end
				if sottoarr=~intpos then		
					restante="#{$'}"
					posizioniarr=restante.split("..")
					posiniziogene=posizioniarr[0]
					if posizioniarr[1]=~strn then	
						posfinegene="#{$`}"
					end
				end	
				posiniziogene=posiniziogene.to_i
				posfinegene=posfinegene.to_i
				posfineplasmo=posfineplasmo.to_i
				posinizioplasmo=posinizioplasmo.to_i
				if ((posiniziogene-posinizioplasmo)>=3) or ((posinizioplasmo-posiniziogene)>=3) then
					if (posfinegene-posiniziogene)==(posfineplasmo-posinizioplasmo) then
						casiparticolari=casiparticolari+k+ " le posizini d'inizio differiscono per almeno 3 basi ma hanno stessa lunghezza\n"
					else
						casiparticolari=casiparticolari+k+" le posizini d'inizio differiscono per almeno 3 basi e non hanno stessa lunghezza\n"
					end
				elsif ((posfineplasmo-posfinegene)>=3) or ((posfinegene-posfineplasmo)>=3) then
					if (posfinegene-posiniziogene)==(posfineplasmo-posinizioplasmo) then
						casiparticolari=casiparticolari+k+ " le posizini di fine differiscono per almeno 3 basi ma hanno stessa lunghezza\n"
					else
						casiparticolari=casiparticolari+k+" le posizini di fine differiscono per almeno 3 basi e non hanno stessa lunghezza\n"
					end
				elsif (posfinegene-posiniziogene)!=(posfineplasmo-posinizioplasmo) then
					casiparticolari=casiparticolari+k+" non hanno stessa lunghezza ma differiscono per meno di 3 basi\n"
				end
			end		
		elsif trasprob!=v[2] then
			casiparticolari=casiparticolari+k+" non trovata isoforma in genedb\n"
		end	
	end
	g.each do |k,v|
		if p.has_key?(k)==false then
			casiparticolari=casiparticolari+k+" non trovato in plasmodb\n"
		end	
	end
end
			
File.open ARGV[1], "w" do |output|
    output << casiparticolari	
end

File.open ARGV[2], "w" do |output|
    output << arraygene
end

#salvo in json i dati dei 14 cromosomi per successiva analisi
dati=arraygene.to_json
File.open('marsdati.txt', 'w') {|f| f.write(dati) }