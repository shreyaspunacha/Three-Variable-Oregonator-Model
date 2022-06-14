run:main.o PhaseField_InitialCondition.o PhaseFieldGeneration.o ThreeVariableOregonatorModel.o 
	gcc -g -o run main.o PhaseField_InitialCondition.o PhaseFieldGeneration.o ThreeVariableOregonatorModel.o -lm 
main.o:main.c
	gcc -g -o main.o -c main.c -lm
ThreeVariableOregonatorModel.o:ThreeVariableOregonatorModel.c
	gcc -g -o ThreeVariableOregonatorModel.o -c ThreeVariableOregonatorModel.c -lm
PhaseFieldGeneration.o:PhaseFieldGeneration.c
	gcc -g -o PhaseFieldGeneration.o -c PhaseFieldGeneration.c
PhaseField_InitialCondition.o:PhaseField_InitialCondition.c
	gcc -g -o PhaseField_InitialCondition.o -c PhaseField_InitialCondition.c 
clean:
	rm main.o PhaseField_InitialCondition.o PhaseFieldGeneration.o ThreeVariableOregonatorModel.o run
