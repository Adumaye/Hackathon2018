# Compilateur Utilisé
CC = mpic++

# Options en mode optimisé - La variable DEBUG est définie comme fausse
OPTIM_FLAG = -O3 -DNDEBUG -w -I Eigen/Eigen -std=c++11

# Options en mode debug - La variable DEBUG est définie comme vraie
DEBUG_FLAG = -g3 -DDEBUG  -I Eigen/Eigen -ltiff -lm -lpthread -std=c++11 -w

# Librairies à linker (création executable)
LIB = -ltiff -lm -lpthread

# On choisit comment on compile
CXX_FLAGS = $(OPTIM_FLAG)

# Le nom de l'exécutable
PROGFilter = mainFilter
PROGSegmentation = mainSegmentation

# Les fichiers source à compiler
SRC = LevelSet.cpp InitMask.cpp ChanVeseSchemes.cpp Image.cpp Util.cpp LevelSet_v.cpp
SRCMainFilter = mainFilter.cc
SRCMainSegmen = mainSegmentation.cc
SRCCompilFilter = mainFilter.o LevelSet.o InitMask.o ChanVeseSchemes.o Image.o Util.o
SRCCompilSegmen = mainSegmentation.o LevelSet.o InitMask.o ChanVeseSchemes.o Image.o Util.o LevelSet_v.o

# La commande complète : compile seulement si un fichier a été modifié
$(PROGSegmentation) : $(SRC) $(SRCMainSegmen)
	$(CC) -c $(SRC) $(SRCMainSegmen) $(CXX_FLAGS)
	$(CC) -o $(PROGSegmentation) $(SRCCompilSegmen) $(LIB)

$(PROGFilter) : $(SRC) $(SRCMainFilter)
	$(CC) -c $(SRC) $(SRCMainFilter) $(CXX_FLAGS)
	$(CC) -o $(PROGFilter) $(SRCCompilFilter) $(LIB)

# Évite de devoir connaitre le nom de l'exécutable
segment : $(PROGSegmentation)
filter : $(PROGFilter)
all : $(PROGFilter) $(PROGSegmentation)

# Supprime l'exécutable, les fichiers binaires (.o) et les fichiers
# temporaires de sauvegarde (~)
# Supprime aussi les fichiers fileint1_aft_preprocess_filtered.tiff fileint1_aft_preprocess_filtered_distance_mask.vtk et fileint1_aft_preprocess_filtered_with_contour.tiff
# (ce qui évite d'avoir à les supprimer soit même lorsque l'on veut tester le code)
clean :
	rm -f *.o *~ $(PROGFilter) *~ $(PROGSegmentation)
	rm -f ../Images/fileint1_aft_preprocess_*
