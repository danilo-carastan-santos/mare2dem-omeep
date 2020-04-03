# mare2dem-omeep
MARE2DEM source files customized for the OMEEP project

Para compilar, é nesessário apenas executar o comando make CLUSTER=omeep
INTEL_PATH=<intel_path>, na pasta raiz do repositório, onde <intel_path> e
o caminho para a pasta intel, contendo os compiladores Intel.
Na pagina de download do MARE2DEM (https://mare2dem.ucsd.edu/?page_id=108) 
existe um conjunto de dados de teste, cujo arquivo de entrada para inversão CSEM 
consiste em leituras referentes à 25 transmissores, 25 receptores e 2 frequências. 
Esse arquivo de entradaencontra-se no diretório invesrion_CSEM. 
Para executar o texte, basta executar o comando mpirun -n 2 <path>/MARE2DEM Demo.0.resistivity dentro
do diretório invesrion_CSEM.
