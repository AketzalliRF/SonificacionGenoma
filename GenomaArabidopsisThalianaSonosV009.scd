///LECTURA DE LOS ARCHIVOS DE TEXTO QUE SE VAN A CARGAR///
(
~dirRoot = "/Users/Aketzalli/Dropbox/Aketzalli/";
//~dirRoot = "/Users/hugosg/Dropbox/UNAM/Aketzalli/";
~dirData = ~dirRoot ++ "data/";
~dirCodigo = ~dirRoot ++ "codigo/";
)
///TERMINA LECTURA DE LOS ARCHIVOS///

(

///INICIA CARGA DE DICCIONARIOS NECESARIOS PARA TRADUCCIÓN ///
~codonAmino = Object.readArchive(~dirCodigo ++ "codonaminoacido.dicc"); // diccionario de la traducción: Codon - aa
~sesgos = Object.readArchive(~dirCodigo ++ "codon_sesgo.dicc");  //diccionario del sesgo de codones
~sesgosTmp = Object.readArchive(~dirCodigo ++ "codon_sesgo.dicc");
~sesgosC = ~sesgosTmp.values.sort.reverse.collect{arg x;
	var objeto = ~sesgosTmp.findKeyForValue(x);
	~sesgosTmp.removeAt(objeto);
	objeto;
}; //ordenar de mayor a menor los codones según la aparición
~sesgosN = ~sesgos.values.reverse.sort.collect{arg x; x}.normalizeSum; //ordenado de mayor a menor
~ocurrencias = Object.readArchive(~dirCodigo ++ "aa_ocurrencia.dicc"); // diccionario de las ocurrencia de aminoacidos de las prot
///TERMINA CARGA DE DICCIONARIOS NECESARIOS PARA TRADUCCIÓN ////////

~genesCompletos = Dictionary(); //los genes completos se refiere a los genes con regiones codificadoras y reguladoras
~genesCDS = Dictionary(); //genes con la región codificadora
~cromosomas = []; //los cromosomas de Arabidopsis thaliana

~leerGenesCompletos = {
	var carpeta = ~dirData ++ "GenesCompletos/";
	carpeta.postln;
	(carpeta ++ "*.txt").pathMatch.do{arg archivo;
		var file, name, contenido;
		name = archivo.split($/).last.split($.)[0].postln;
		file = FileReader.read(archivo);
		contenido = file[1..].flat;
		~genesCompletos.put(name.toUpper.asSymbol, Dictionary().put(\bases, contenido))
	};
}.value;

~leerGenesCDS = {
	var carpeta = ~dirData ++ "GenesCDS/";
	carpeta.postln;
	(carpeta ++ "*.txt").pathMatch.do{arg archivo;
		var file, name, contenido;
		name = archivo.split($/).last.split($.)[0].postln;
		file = FileReader.read(archivo);
		contenido = file[1..].flat;
		~genesCDS.put(name.toUpper.asSymbol, Dictionary().put(\bases, contenido))
	};
}.value;
 ///Codones///
~generarCodones = {
	~genesCDS.do{arg item;
		item.at(\bases).clump(3);
		item.put(\codones, item.at(\bases).clump(3).collect{arg each;
			(each[0] ++ each[1] ++ each[2]).asSymbol;
		});
	};
}.value;

///Sesgo del codon///
~generarSesgo = {
	~genesCDS.do{arg item;
		item.at(\codones);
		item.put(\sesgo, item.at(\codones).collect{arg x;
			~sesgos.at(x.asSymbol);
		});
	};
}.value;

///Traducción del codigo genetico///
~generarAminoacidos = {
	~genesCDS.do{arg item;
		item.at(\codones).postln;
		item.put(\aminoacidos, item.at(\codones).collect{arg x;
			~codonAmino.at(x.asSymbol);
		});
	};
}.value;

///Ocurrencias de los aminoacidos en las prteinas////
~generarOcurrencias = {
	~genesCDS.do{arg item;
		item.at(\aminoacidos);
		item.put(\ocurrencias, item.at(\aminoacidos).collect{arg x;
			~ocurrencias.at(x.asSymbol);
		});
	};
}.value;

///Cromosomas///
~leerCromosomas = {
	var carpeta = ~dirData ++ "genoma/";
	carpeta.postln;
	5.do{arg crom;
		var file;
		file = FileReader.read(carpeta ++ (crom + 1) ++ ".txt").flat;
		file.postln;
		~cromosomas = ~cromosomas.add(file);
	}
}.value;
)

///INICIO EJEMPLOS DE USO ///
~codonesDia = ~genesCDS.at(\DIA).at(\bases)
~codonesDia = ~genesCDS.at(\FLC).at(\codones)
~codonesDia = ~genesCDS.at(\FLC).at(\sesgo)
~codonesDia = ~genesCDS.at(\FLC).at(\aminoacidos).size
~codonesDia = ~genesCDS.at(\FLC).at(\ocurrencias).normalize
///FINAL EJEMPLOS DE USO ////////

///INICIO SONIFICACION DE GENOMA COMPLETO///
SynthDef(\basenitro, {arg freq,bati, vol, pan=1, dur;
	var signal, env;
	env  = EnvGen.kr(Env.sine(dur), doneAction: 2);
	signal = SinOsc.ar([freq,freq+bati]).sum;
	Out.ar(0, Pan2.ar(signal * vol, pan, env));
}).add;


//5 cromosomas en A. thaliana
(
Ppar([
	Pbind(
	\instrument, \basenitro,
	\bases, Pseq(~cromosomas[0], inf).trace, //Cromosoma 1
	\dur, Pfunc{|ev| switch(ev[\bases],$A,{1},$G,{6},$T,{4},$C,{2},$N,{0.5},{1})}, //posición del C=O que forma puente de hidrógeno
	\delta, Pfunc{|ev| switch(ev[\bases],$G,{2},$C,{1},$T,{1},$A,{2},$N,{0.5},{1})}, //la cantidad de anillos de la estructura química
	\pan, Pfunc{|ev| switch(ev[\bases],$A,{1},$G,{0.5},$T,{-0.5},$C,{-1},$N,{0},{0})},
	\vol, Pfunc{|ev| switch(ev[\bases],$A,{0.2},$G,{0.3},$T,{0.2},$C,{0.3},$N,{0.01},{0.01}) * 0.15}, //cuantos puentes de hidrógeno forman
	\freq, Pfunc{|ev| switch(ev[\bases],$A,{[1,6]},$G,{[1,2]},$T,{3},$C,{[3,4]},$N,{0},{0})*60}, //posición de los N que forman los puente de hidrógeno
	\bati, Pfunc{|ev| switch(ev[\bases],$A,{([1,6]).choose},$G,{([1,2,6]).choose},$T,{([3,4]).choose},$C,{([2,3,4]).choose},$N,{0},{0})*1.5}, //posición de todos los puentes de hidrógeno
),

	Pbind(
	\instrument, \basenitro,
	\bases, Pseq(~cromosomas[1], inf).trace, //Cromosoma 2
	\dur, Pfunc{|ev| switch(ev[\bases],$A,{1},$G,{6},$T,{4},$C,{2},$N,{0.5},{1})},
	\delta, Pfunc{|ev| switch(ev[\bases],$G,{2},$C,{1},$T,{1},$A,{2},$N,{0.5},{1})},
	\pan, Pfunc{|ev| switch(ev[\bases],$A,{1},$G,{0.5},$T,{-0.5},$C,{-1},$N,{0},{0})},
	\vol, Pfunc{|ev| switch(ev[\bases],$A,{0.2},$G,{0.3},$T,{0.2},$C,{0.3},$N,{0.01},{0.01}) * 0.1},
	\freq, Pfunc{|ev| switch(ev[\bases],$A,{[1,6]},$G,{[1,2]},$T,{3},$C,{[3,4]},$N,{0},{0})*70},
	\bati, Pfunc{|ev| switch(ev[\bases],$A,{([1,6]).choose},$G,{([1,2,6]).choose},$T,{([3,4]).choose},$C,{([2,3,4]).choose},$N,{0},{0})*1.5},
),

	Pbind(
	\instrument, \basenitro,
	\bases, Pseq(~cromosomas[2], inf).trace, //Cromosoma 3
	\dur, Pfunc{|ev| switch(ev[\bases],$A,{1},$G,{6},$T,{4},$C,{2},$N,{0.5},{1})},
	\delta, Pfunc{|ev| switch(ev[\bases],$G,{2},$C,{1},$T,{1},$A,{2},$N,{0.5},{1})},
	\pan, Pfunc{|ev| switch(ev[\bases],$A,{1},$G,{0.5},$T,{-0.5},$C,{-1},$N,{0},{0})},
	\vol, Pfunc{|ev| switch(ev[\bases],$A,{0.2},$G,{0.3},$T,{0.2},$C,{0.3},$N,{0.01},{0.01}) * 0.1},
	\freq, Pfunc{|ev| switch(ev[\bases],$A,{[1,6]},$G,{[1,2]},$T,{3},$C,{[3,4]},$N,{0},{0})*80},
	\bati, Pfunc{|ev| switch(ev[\bases],$A,{([1,6]).choose},$G,{([1,2,6]).choose},$T,{([3,4]).choose},$C,{([2,3,4]).choose},$N,{0},{0})*1.5},
),
	Pbind(
	\instrument, \basenitro,
	\bases, Pseq(~cromosomas[3], inf).trace, //Cromosoma 4
	\dur, Pfunc{|ev| switch(ev[\bases],$A,{1},$G,{6},$T,{4},$C,{2},$N,{0.5},{1})},
	\delta, Pfunc{|ev| switch(ev[\bases],$G,{2},$C,{1},$T,{1},$A,{2},$N,{0.5},{1})},
	\pan, Pfunc{|ev| switch(ev[\bases],$A,{1},$G,{0.5},$T,{-0.5},$C,{-1},$N,{0},{0})},
	\vol, Pfunc{|ev| switch(ev[\bases],$A,{0.2},$G,{0.3},$T,{0.2},$C,{0.3},$N,{0.01},{0.01}) * 0.1},
	\freq, Pfunc{|ev| switch(ev[\bases],$A,{[1,6]},$G,{[1,2]},$T,{3},$C,{[3,4]},$N,{0},{0})*90},
	\bati, Pfunc{|ev| switch(ev[\bases],$A,{([1,6]).choose},$G,{([1,2,6]).choose},$T,{([3,4]).choose},$C,{([2,3,4]).choose},$N,{0},{0})*1.5},
),
	Pbind(
	\instrument, \basenitro,
	\bases, Pseq(~cromosomas[4], inf).trace, //Cromosoma 5
	\dur, Pfunc{|ev| switch(ev[\bases],$A,{1},$G,{6},$T,{4},$C,{2},$N,{0.5},{1})},
	\delta, Pfunc{|ev| switch(ev[\bases],$G,{2},$C,{1},$T,{1},$A,{2},$N,{0.5},{1})},
	\pan, Pfunc{|ev| switch(ev[\bases],$A,{1},$G,{0.5},$T,{-0.5},$C,{-1},$N,{0},{0})},
	\vol, Pfunc{|ev| switch(ev[\bases],$A,{0.2},$G,{0.3},$T,{0.2},$C,{0.3},$N,{0.01},{0.01}) * 0.1},
	\freq, Pfunc{|ev| switch(ev[\bases],$A,{[1,6]},$G,{[1,2]},$T,{3},$C,{[3,4]},$N,{0},{0})*100},
	\bati, Pfunc{|ev| switch(ev[\bases],$A,{([1,6]).choose},$G,{([1,2,6]).choose},$T,{([3,4]).choose},$C,{([2,3,4]).choose},$N,{0},{0})*1.5},
),
],inf).play
)
///FINAL SONIFICACION DE GENOMA COMPLETO///

///PRUEBA DE BOCINAS///
{Mix.new (SinOsc.ar([600]))}.scope;
play({Pan2.ar( SinOsc.ar(XLine.kr(20,5000,10),0,0.5),-1 )}); //-1 izq, 1 derecha.
///FIN DE PRUEBA DE BOCINAS///


