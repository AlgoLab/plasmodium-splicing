# encoding: utf-8

require 'fileutils'
require 'json'

#funzione che calcola la posizione di inizio del primo esone tra le varie isoforme 
def minimo (array)
	min=4000000
	array.each do |a|
		if a.to_i<min then
			min=a.to_i
		end
	end
	return min
end

#funzione che calcola la posizione di fine per l'ultimo esone tra le varie isoforme 
def massimo (array)
	max=0
	array.each do |a|
		if a.to_i>max then
			max=a.to_i
		end
	end
	return max
end

#ogni volta che viene rieseguito il programma viene cancellata la vecchia gerarchia di cartelle contenenti i dati
1.upto(14) do |n|
    if Dir.exists?('crom'+n.to_s) then
	  FileUtils.remove_dir('crom'+n.to_s)
    end
end

#dopo aver analizzato i due db viene recuperata la struttura dati che era stata creata
app=File.read("marsdati.txt")
#e poi salvata
hashcompleta = JSON::parse(app)

#variabili utilizzate
punto=/\.1/
due=/\.2/
sq=/^\s+/
puntotre=/\.3/
parentesiap=/\(/
parentesich=/\)/
posdue=[]
posuno=[]
postre=[]
arraytre=[]
estrai=""
nomeprimicrom="Pf3D7_0"
nomecrom="Pf3D7_"
finenomecrom="_v3.embl"
arraysequenze=Array.new(14)
strmeno=Array.new(14)
hashstampa = Hash.new
strprovv=""
ric=[]
arr=[]
val=[]

#per ogni cromosoma
1.upto(14) do |n|
	#viene creata la cartella corrispondente
	Dir.mkdir('crom'+n.to_s)
	if n<10 then
		nomefile=nomeprimicrom+n.to_s+finenomecrom
	else
		nomefile=nomecrom+n.to_s+finenomecrom
	end   
	#viene poi cercato l'embl corrispondente
	File.open(nomefile, 'r') do |leggi|
 		leggi.each_line do |linea|
            linea.chomp!
			if linea=~sq then
				#ed infine estratta la sequenza
				strprovv=strprovv+linea
			end
		end			
	end
	#vengono eliminati i numeri finali di ogni riga e gli spazi che ci sono ogni 10 basi
	strprovv=strprovv.delete " 1234567890/"
	#una volta ripulita viene salvata in un array che poi conterra tutte le sequenze dei 14 cromosomi
	arraysequenze[n-1]=strprovv
	#viene invertita 
	nuova=strprovv.reverse!    
	#e complementata
    nuova.gsub!(/[acgt]/, 'a' => 't', 'c' => 'g', 't' => 'a', 'g' => 'c')
	#cosi viene gia generato l'array con le 14 sequenze soggette a reverse&complement per gli strand -
    strmeno[n-1]=nuova
	strprovv=""
end
sottostr=""

#per generare i file
1.upto(14) do |n|
	chiavi=[]
	nuovonome=""
	dastamp=""
	#seleziono i dati inerenti a quel cromosoma
	hashgenomiche=Hash.new
	hashmanc=Hash.new
	arr=hashcompleta[n-1]
	#accedo alla cartella del cromosoma
    if n!=1 then
        Dir.chdir("..")
        Dir.chdir("crom#{n}")
    else
    	Dir.chdir("crom#{n}")
    end
	genom=""
	#creo un array contenente tutte le chiavi
	arr.each do |k, val|
		chiavi.push(k)
	end
	#e lo ordino
	chiavi.sort!
	#per ognuna delle chiavi 
	chiavi.each do |k|
		val=arr[k]
		estrai=val[4] 
		#estraggo il segno della isoforma (ho in precedenza controllato che le diverse isoforme avessero stesso segno)
		if estrai=~parentesiap then
			rest="#{$'}" 
			if rest=~parentesich then
				segno="#{$`}"
			end
		end 
		sottosq=""
		inizio=0
		fine=0	
		ant=""
		#rimuovo il . dal nome univoco per pintron
		nuovo=k.split(".")
		if k=~punto then
			if segno=="+" then
				sottostr=arraysequenze[n-1]
			else
				sottostr=strmeno[n-1]
			end
			Dir.mkdir(k)
			nuovonome=nuovo[0]+"_"+nuovo[1]
			posizioni=val[5]        
			percorso=k+"/"+nuovonome+"genom.txt"        
			if posizioni[0].instance_of?(String) then
				lung=posizioni.last.to_i-posizioni.first.to_i+2001				
				inizio=posizioni.first.to_i-1001
				fine=posizioni.last.to_i+999
				sottosq=sottostr[inizio..fine]
			else
				lung=posizioni.last.last.to_i-posizioni.first.first.to_i+2001				
				inizio=posizioni.first.first.to_i-1001
				fine=posizioni.last.last.to_i+999
				sottosq=sottostr[inizio..fine]
			end        
			dastamp=">chr1:#{inizio+1}:#{fine+1}:+1\n"
			genom=dastamp+sottosq
			hashgenomiche[k]=genom            
			Dir.chdir(k)
			File.open(nuovonome+"genom.txt", 'w') do |output|
				output << genom
			end
			Dir.chdir("..")
		else
			ant=k[0..12]
			antenato=ant+".1"
			nuovonome=nuovo[0]
			nuovonome=nuovo[0]+"_"+nuovo[1]
			posizioni=val[5]        
			percorso=k+"/"+nuovonome+"genom.txt"
			contant=arr[antenato]
			estrai=contant[4]
			if estrai=~parentesiap then
				rest="#{$'}" 
				if rest=~parentesich then
					segnoantenato="#{$`}"
				end
			end
			hashmanc[k]=[antenato,percorso,ant,segno,segnoantenato]
		end
		
	end
    i=0
	f=0
	danonrifare=[]
	#per ogni gene che presenta più isoforme
	hashmanc.each do |k, val|
		percorso=val[1]
		antenato=val[0]
		nn=val[2]
		s=val[3]
		if k=~puntotre then
			#prendo in considerazione le altre due isoforme
			v=arr[k]
			postre=v[5]
			#rimuovo la cartella che era stata creata per il .1 in quanto non è piu valida 
			FileUtils.remove_dir(antenato)
			#e ne creo una nuova generica che conterra le diverse isoforme
			Dir.mkdir(nn)
			andue=nn+".2"
			arraytre.push(andue)
			valdue=arr[andue]
			danonrifare.push(valdue)
			valant=arr[antenato]
			ric.push(antenato)
			posdue=valdue[5]
			posuno=valant[5]
			postre=postre.flatten
			posdue=posdue.flatten
			posuno=posuno.flatten
			posuno.concat(posdue)
			posuno.concat(postre)
			i=minimo(posuno)-1000
			f=massimo(posuno)+1000
			isodue=hashmanc[andue]
			isouno=val[4]
			sdue=isodue[3]
			# if (s!=sdue) or (s!=isouno) then
				# puts "errore segno isoforme"
				# puts k
			# end
			if sdue=="-" then
				dastamp=">chr1:#{i}:#{f}:-1\n"
				gencompl=dastamp+sottostr[i-1..f-1]
			else
				dastamp=">chr1:#{i}:#{f}:+1\n"
				gencompl=dastamp+sottostr[i-1..f-1]
			end
			Dir.chdir(nn)
			File.open(nn+"genom.txt", 'w') do |output|
				output << gencompl
			end
			Dir.chdir("..")
		end
		if (k=~due) and (danonrifare.include?(k)==false) and (hashmanc.include?(nn+".3")==false) then
			FileUtils.remove_dir(antenato)
			ric.push(antenato)
			Dir.mkdir(nn)
			v=arr[k]
			posdue=v[5]
			valant=arr[antenato]
			posuno==valant[5]
			posdue=posdue.flatten(1)
			posuno=posuno.flatten
			posuno.concat(posdue)
			i=minimo(posuno)-1000
			f=massimo(posuno)+1000
			isouno=val[4]			
			# if s!=isouno then
				# puts "errore segno isoforme"
				# puts k
			# end
			if sdue=="-" then
				dastamp=">chr1:#{i}:#{f}:-1\n"
				gencompl=dastamp+sottostr[i-1..f-1]
			else
				dastamp=">chr1:#{i}:#{f}:+1\n"
				gencompl=dastamp+sottostr[i-1..f-1]
			end
			Dir.chdir(nn)
			File.open(nn+"genom.txt", 'w') do |output|
				output << gencompl
			end
			Dir.chdir("..")
		end
		if sdue=="-" then
				dastamp=">chr1:#{i}:#{f}:-1\n"
				gencompl=dastamp+sottostr[i-1..f-1]
		else
				dastamp=">chr1:#{i}:#{f}:+1\n"
				gencompl=dastamp+sottostr[i-1..f-1]
		end
		hashstampa[nn]=[nuovonome, gencompl]		
	end
	
end

Dir.chdir("..")

#per creare i file con cds
1.upto(14) do |n|
	chiavi=[]
	if n!=1 then
        Dir.chdir("..")
        Dir.chdir("crom#{n}")
    else
    	Dir.chdir("crom#{n}")
    end
	arr=hashcompleta[n-1]
	arr.each do |k, val|
		chiavi.push(k)
	end
	chiavi.sort!
	chiavi.each do |k|
		analizzati=[]
		cds=""
		val=arr[k]
		estrai=val[4]
		if estrai=~parentesiap then
			rest="#{$'}" 
			if rest=~parentesich then
				segno="#{$`}"
			end
		end 
		if segno=="+" then
			sottostr=arraysequenze[n-1]
		else
			sottostr=strmeno[n-1]
		end
		nuovo=k.split(".")
		nuovonome=nuovo[0]+"_"+nuovo[1]	
		stampa=">/gb="+nuovonome+" /clone_end=3'\n"
		if (k=~punto) and (ric.include?(k)==false) then
			Dir.chdir(k)
			posizioni=val[5]
			posizioni.each do |sottoarr|
				if sottoarr.instance_of?(String) then
					cds=cds+sottostr[sottoarr.to_i..sottoarr.to_i]
				else
					cds=cds+sottostr[sottoarr.first.to_i..sottoarr.last.to_i]
				end
			end
			stampa=stampa+cds
			File.open(k+".txt", 'w') do |output|
				output << stampa
			end
		else
			ant=k[0..12]
			Dir.chdir(ant)
			if (k=~due) and (arraytre.include?(k)==false) then
				cds=">/gb="+ant+"_1 /clone_end=3'\n"
				posizionidue=val[5]
				p=ant+".1"
				valuno=arr[p]
				posizioniuno=valuno[5]
				posizioniuno.each do |sottoarr|
					cds=cds+sottostr[sottoarr.first.to_i..sottoarr.last.to_i]
				end
				cds=cds+"\n>/gb="+nuovonome+" /clone_end=3'\n"
				posizionidue.each do |sottoarr|
					cds=cds+sottostr[sottoarr.first.to_i..sottoarr.last.to_i]
				end
				
			else
				cds=">/gb="+ant+"_1 /clone_end=3'\n"
				posizionitre=val[5]
				p=ant+".1"
				d=ant+".2"
				valuno=arr[p]
				valdue=arr[d] 
				posizionidue=valdue[5]
				posizioniuno=valuno[5]
				posizioniuno.each do |sottoarr|
					cds=cds+sottostr[sottoarr.first.to_i..sottoarr.last.to_i]
				end
				cds=cds+"\n>/gb="+ant+"_2 /clone_end=3'\n"
				posizionidue.each do |sottoarr|
					cds=cds+sottostr[sottoarr.first.to_i..sottoarr.last.to_i]
				end
				cds=cds+"\n>/gb="+nuovonome+" /clone_end=3'\n"
				posizionitre.each do |sottoarr|
					cds=cds+sottostr[sottoarr.first.to_i..sottoarr.last.to_i]
				end
				
			end
			File.open(ant+".txt", 'w') do |output|
				output << cds
			end
		end
		Dir.chdir("..")
	end
end
