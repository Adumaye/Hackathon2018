# Compilateur Utilisé
CC = mpic++

# Options en mode optimisé - La variable DEBUG est définie comme fausse
OPTIM_FLAG = -O3 -DNDEBUG -w -I Eigen/Eigen -std=c++11

# Options en mode debug - La variable DEBUG est définie comme vraie
DEBUG_FLAG = -g -DDEBUG  -I Eigen/Eigen -ltiff -lm -lpthread -std=c++11 -w

# Librairies à linker (création executable)
LIB = -ltiff -lm -lpthread

# On choisit comment on compile
CXX_FLAGS = $(OPTIM_FLAG)

PLAFRIM_FLAG = -DNDEBUG -ta=tesla:managed -I ../eigen/Eigen -Minfo=accel -acc -O3 -w -std=c++11

# Le nom de l'exécutable
PROGFilter = mainFilter
PROGSegmentation = mainSegmentation
PROGFilterPlaf = mainFilterPlaf
PROGSegmentationPlaf = mainSegmentationPlaf
# Les fichiers source à compiler
SRC = InitMask.cpp ChanVeseSchemes.cpp Image.cpp Util.cpp LevelSet_v.cpp
SRCMainFilter = mainFilter.cc
SRCMainSegmen = mainSegmentation.cc
SRCCompilFilter = mainFilter.o InitMask.o ChanVeseSchemes.o Image.o Util.o LevelSet_v.o
SRCCompilSegmen = mainSegmentation.o InitMask.o ChanVeseSchemes.o Image.o Util.o LevelSet_v.o

# La commande complète : compile seulement si un fichier a été modifié
$(PROGSegmentation) : $(SRC) $(SRCMainSegmen)
	$(CC) -c $(SRC) $(SRCMainSegmen) $(CXX_FLAGS)
	$(CC) -o $(PROGSegmentation) $(SRCCompilSegmen) $(LIB)

$(PROGFilter) : $(SRC) $(SRCMainFilter)
	$(CC) -c $(SRC) $(SRCMainFilter) $(CXX_FLAGS)
	$(CC) -o $(PROGFilter) $(SRCCompilFilter) $(LIB)

#Pour PlafRIM
$(PROGSegmentationPlaf) : $(SRC) $(SRCMainSegmen)
	$(CC) $(PLAFRIM_FLAG) -c $(SRC) $(SRCMainSegmen)
	$(CC) -o $(PROGSegmentation) $(SRCCompilSegmen) $(LIB) -ta=tesla:managed

$(PROGFilterPlaf) : $(SRC) $(SRCMainFilter)
	$(CC) $(PLAFRIM_FLAG) -c $(SRC) $(SRCMainFilter)
	$(CC) -o $(PROGFilter) $(SRCCompilFilter) $(LIB) -ta=tesla:managed

# Évite de devoir connaitre le nom de l'exécutable
all : $(PROGFilter)     $(PROGSegmentation)
filter : $(PROGFilter)
segment : $(PROGSegmentation)

plafrim : $(PROGFilterPlaf)     $(PROGSegmentationPlaf)

# Supprime l'exécutable, les fichiers binaires (.o) et les fichiers
# temporaires de sauvegarde (~)
clean :
	rm -f *.o *~ $(PROGFilter) *~ $(PROGSegmentation)

cleanall :
	rm -f *.o *~ $(PROGFilter) *~ $(PROGSegmentation)
	#rm -rf ../Images/fileint1_aft_preprocess_filtered.tiff
	rm -rf ../Images/fileint1_aft_preprocess_filtered_distance_mask.vtk
