#!/bin/csh
#		Prepare additional codon usage tables

$BLIMPS_DIR/bin/coduse CUTG/gbpln.spsum "Arabidopsis thaliana" arab.codon.use
$BLIMPS_DIR/bin/coduse CUTG/gbpln.spsum "Brassica oleracea" broc.codon.use
$BLIMPS_DIR/bin/coduse CUTG/gbpln.spsum "Brassica napus" napus.codon.use
$BLIMPS_DIR/bin/coduse CUTG/gbpln.spsum "Oryza sativa" rice.codon.use
$BLIMPS_DIR/bin/coduse CUTG/gbpln.spsum "Nicotiana tabacum" tobacco.codon.use
$BLIMPS_DIR/bin/coduse CUTG/gbpln.spsum "Zea mays" zea.codon.use
$BLIMPS_DIR/bin/coduse CUTG/gbpln.spsum "Hordeum vulgare" barley.codon.use
$BLIMPS_DIR/bin/coduse CUTG/gbpln.spsum "Lycopersicon esculentum" tomato.codon.use
$BLIMPS_DIR/bin/coduse CUTG/gbpln.spsum "Neurospora crassa" neur.codon.use
$BLIMPS_DIR/bin/coduse CUTG/gbpln.spsum "Chlamydomonas reinhardtii" chlamy.codon.use
$BLIMPS_DIR/bin/coduse CUTG/gbpln.spsum "Schizosaccharomyces pombe" pombe.codon.use
$BLIMPS_DIR/bin/coduse CUTG/gbpln.spsum "Saccharomyces cerevisiae" yeast.codon.use

$BLIMPS_DIR/bin/coduse CUTG/gbbct.spsum "Escherichia coli" ecoli.codon.use
$BLIMPS_DIR/bin/coduse CUTG/gbbct.spsum "Streptomyces coelicolor" strep.codon.use
$BLIMPS_DIR/bin/coduse CUTG/gbbct.spsum "Lactobacillus casei" lacto.codon.use
$BLIMPS_DIR/bin/coduse CUTG/gbbct.spsum "Staphylococcus aureus" staph.codon.use
$BLIMPS_DIR/bin/coduse CUTG/gbbct.spsum "Bartonella bacilliformis" bartb.codon.use
$BLIMPS_DIR/bin/coduse CUTG/gbbct.spsum "Bartonella henselae" barth.codon.use
$BLIMPS_DIR/bin/coduse CUTG/gbbct.spsum "Ehrlichia chaffeensis" ehrl.codon.use
$BLIMPS_DIR/bin/coduse CUTG/gbbct.spsum "Bacillus subtilis" bacsub.codon.use
$BLIMPS_DIR/bin/coduse CUTG/gbbct.spsum "Bacillus thuringiensis" bact.codon.use
$BLIMPS_DIR/bin/coduse CUTG/gbbct.spsum "Agrobacterium tumefaciens" agro.codon.use
$BLIMPS_DIR/bin/coduse CUTG/gbbct.spsum "Chlamydia trachomatis" chlam.codon.use

$BLIMPS_DIR/bin/coduse CUTG/gbinv.spsum "Giardia intestinalis" giardia.codon.use
$BLIMPS_DIR/bin/coduse CUTG/gbinv.spsum "Tetrahymena thermophila" tetra.codon.use
$BLIMPS_DIR/bin/coduse CUTG/gbinv.spsum "Trypanosoma brucei" tryp.codon.use
$BLIMPS_DIR/bin/coduse CUTG/gbinv.spsum "Leishmania tarentolae" leish.codon.use
$BLIMPS_DIR/bin/coduse CUTG/gbinv.spsum "Plasmodium falciparum" plasmo.codon.use
$BLIMPS_DIR/bin/coduse CUTG/gbinv.spsum "Schistosoma mansoni" schis.codon.use
$BLIMPS_DIR/bin/coduse CUTG/gbinv.spsum "Bombyx mori" bombyx.codon.use
$BLIMPS_DIR/bin/coduse CUTG/gbinv.spsum "Caenorhabditis elegans" worm.codon.use
$BLIMPS_DIR/bin/coduse CUTG/gbinv.spsum "Drosophila melanogaster" fly.codon.use
$BLIMPS_DIR/bin/coduse CUTG/gbinv.spsum "Hydra vulgaris" hydra.codon.use
$BLIMPS_DIR/bin/coduse CUTG/gbinv.spsum "Brugia malayi" brugia.codon.use
$BLIMPS_DIR/bin/coduse CUTG/gbinv.spsum "Onchocerca volvulus" onch.codon.use
$BLIMPS_DIR/bin/coduse CUTG/gbinv.spsum "Ascaris suum" ascaris.codon.use
$BLIMPS_DIR/bin/coduse CUTG/gbinv.spsum "Trichinella spiralis" trich.codon.use
$BLIMPS_DIR/bin/coduse CUTG/gbinv.spsum "Haemonchus contortus" haem.codon.use
$BLIMPS_DIR/bin/coduse CUTG/gbinv.spsum "Ancylostoma caninum" ancyl.codon.use
$BLIMPS_DIR/bin/coduse CUTG/gbinv.spsum "Lymnaea stagnalis" lymn.codon.use
$BLIMPS_DIR/bin/coduse CUTG/gbinv.spsum "Placopecten magellanicus" placo.codon.use

$BLIMPS_DIR/bin/coduse CUTG/gbvrt.spsum "Danio rerio" zebra.codon.use
$BLIMPS_DIR/bin/coduse CUTG/gbvrt.spsum "Fugu rubripes" fugu.codon.use
$BLIMPS_DIR/bin/coduse CUTG/gbvrt.spsum "Xenopus laevis" frog.codon.use
$BLIMPS_DIR/bin/coduse CUTG/gbvrt.spsum "Ambystoma mexicanum" amby.codon.use
$BLIMPS_DIR/bin/coduse CUTG/gbvrt.spsum "Gallus gallus" chicken.codon.use
$BLIMPS_DIR/bin/coduse CUTG/gbvrt.spsum "Anas platyrhynchos" anas.codon.use
$BLIMPS_DIR/bin/coduse CUTG/gbvrt.spsum "Cairina moschata" duck.codon.use
$BLIMPS_DIR/bin/coduse CUTG/gbvrt.spsum "Anser anser" goose.codon.use

$BLIMPS_DIR/bin/coduse CUTG/gbrod.spsum "Mus musculus" mouse.codon.use
$BLIMPS_DIR/bin/coduse CUTG/gbrod.spsum "Rattus norvegicus" rat.codon.use
$BLIMPS_DIR/bin/coduse CUTG/gbrod.spsum "Cavia porcellus" cavia.codon.use

$BLIMPS_DIR/bin/coduse CUTG/gbmam.spsum "Bos taurus" cow.codon.use
$BLIMPS_DIR/bin/coduse CUTG/gbmam.spsum "Sus scrofa" pig.codon.use
$BLIMPS_DIR/bin/coduse CUTG/gbmam.spsum "Oryctolagus cuniculus" rabbit.codon.use

$BLIMPS_DIR/bin/coduse CUTG/gbpri.spsum "Homo sapiens" human.codon.use

exit
