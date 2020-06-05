#include <iostream>
#include "math_tools.h"
#include "classes.h"
#include "tools.h"
#include "display_tools.h"
#include "sel.h"
#include "assembly.h"

int main(int argc, char *argv[])
{
    char filename[150];
    strcpy(filename,argv[1]);

    vector<Matrix> localKs;
    vector<Vector> localbs;
    Matrix K;
    Vector b;
    Vector T;

    cout << "IMPLEMENTACI"<<char(224)<<"N DEL M"<<char(144)<<"TODO DE LOS ELEMENTOS FINITOS\n"
         << "\t- ECUACIONES DE NAVIER-STOKES\n" << "\t- 2 DIMENSIONES\n"
         << "\t- FUNCIONES DE FORMA LINEALES\n" << "\t- PESOS DE GALERKIN\n"
         << "*********************************************************************************\n\n";

    mesh m;
    leerMallayCondiciones(m,filename);
    cout << "Datos obtenidos correctamente\n********************\n";

    crearSistemasLocales(m,localKs,localbs);
    showKs(localKs); showbs(localbs);
    cout << "******************************\n";

    zeroes(K,3*m.getSize(NODES));
    zeroes(b,3*m.getSize(NODES));
    ensamblaje(m,localKs,localbs,K,b);
    showMatrix(K); showVector(b);
    cout << "******************************\n";

    applyDirichlet(m,K,b);
    showMatrix(K); showVector(b);
    cout << "******************************\n";

    zeroes(T,b.size());
    calculate(K,b,T);

    cout << "La respuesta es: \n";
    showVector(T);

    writeResults(m,T,filename);

    return 0;
}