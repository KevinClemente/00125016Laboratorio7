#include <fstream>
#include "string.h"

void obtenerDatos(istream &file,int nlines,int n,int mode,item* item_list){
    string line;
    file >> line;
    if(nlines==DOUBLELINE) file >> line;

    for(int i=0;i<n;i++){
        switch(mode){
        case INT_FLOAT:
            int e0; float r0;
            file >> e0 >> r0;
            item_list[i].setValues(0,0,0,e0,0,0,r0);
            break;
        case INT_FLOAT_FLOAT:
            int e; float r,rr;
            file >> e >> r >> rr;
            item_list[i].setValues(e,r,rr,0,0,0,0);
            break;
        case INT_INT_INT_INT:
            int e1,e2,e3,e4;
            file >> e1 >> e2 >> e3 >> e4;
            item_list[i].setValues(e1,0,0,e2,e3,e4,0);
            break;
        }
    }
}

void correctConditions(int n,condition *list,int *indices){
    for(int i=0;i<n;i++)
        indices[i] = list[i].getNode1();

    for(int i=0;i<n-1;i++){
        int pivot = list[i].getNode1();
        for(int j=i;j<n;j++)
            //Si la condici�n actual corresponde a un nodo posterior al nodo eliminado por
            //aplicar la condici�n anterior, se debe actualizar su posici�n.
            if(list[j].getNode1()>pivot)
                list[j].setNode1(list[j].getNode1()-1);
    }
}

void addExtension(char *newfilename,char *filename,char *extension){
    int ori_length = strlen(filename);
    int ext_length = strlen(extension);
    int i;
    for(i=0;i<ori_length;i++)
        newfilename[i] = filename[i];
    for(i=0;i<ext_length;i++)
        newfilename[ori_length+i] = extension[i];
    newfilename[ori_length+i] = '\0';
}

void correctIndices(int n,condition* list,int total){
    for(int i=0;i<n;i++)
        list[i].setNode1(list[i].getNode1()+total);
}

void addArray(int* index, int n,condition* list, condition* listf){
    for(int i=0;i<n;i++){
        listf[*index] = list[i];
        (*index)++;
    }
}

void fusionDirichlet(int n1,condition* list1,int n2,condition* list2,int n3,condition* list3,condition* list){
    int index = 0;
    addArray(&index,n1,list1,list);
    addArray(&index,n2,list2,list);
    addArray(&index,n3,list3,list);
}

void leerMallayCondiciones(mesh &m,char *filename){
    char inputfilename[150];
    ifstream file;
    float a,b,c,d,e,f,g;
    int nnodes,neltos,ndirich_u,ndirich_v,ndirich_p;
    condition *dirichlet_u, *dirichlet_v, *dirichlet_p;

    addExtension(inputfilename,filename,".dat");
    file.open(inputfilename);

    file >> a >> b >> c >> d >> e >> f >> g;
    //cout <<u_bar<<nu<<rho<<f_x<<f_y<<"\n";
    file >> nnodes >> neltos >> ndirich_u >> ndirich_v >>ndirich_p;
    //cout <<nnodes<<neltos<<ndirich_u<<ndirich_v<<ndirich_p<<"\n";

    m.setParameters(a,b,c,d,e,f,g);
    m.setSizes(nnodes,neltos,ndirich_u+ndirich_v+ndirich_p);
    m.createData();

    dirichlet_u = new condition[ndirich_u];
    dirichlet_v = new condition[ndirich_v];
    dirichlet_p = new condition[ndirich_p];

    obtenerDatos(file,SINGLELINE,nnodes,INT_FLOAT_FLOAT,m.getNodes());
    obtenerDatos(file,DOUBLELINE,neltos,INT_INT_INT_INT,m.getElements());
    obtenerDatos(file,DOUBLELINE,ndirich_u,INT_FLOAT,dirichlet_u);
    obtenerDatos(file,DOUBLELINE,ndirich_v,INT_FLOAT,dirichlet_v);
    obtenerDatos(file,DOUBLELINE,ndirich_p,INT_FLOAT,dirichlet_p);

    file.close();

    correctIndices(ndirich_v,dirichlet_v,nnodes);
    correctIndices(ndirich_p,dirichlet_p,2*nnodes);

    fusionDirichlet(ndirich_u,dirichlet_u,ndirich_v,dirichlet_v,ndirich_p,dirichlet_p,m.getDirichlet());

    correctConditions(ndirich_u+ndirich_v+ndirich_p,m.getDirichlet(),m.getDirichletIndices());
}

bool findIndex(int v, int s, int *arr){
    for(int i=0;i<s;i++)
        if(arr[i]==v) return true;
    return false;
}

int getIndex(int v, int s, int *arr){
    for(int i=0;i<s;i++)
        if(arr[i]==v) return i;
    return -1;
}

int getIndex(int v, int s, Vector vec){
    for(int i=0;i<s;i++)
        if(vec.at(i)==v) return i;
    return -1;
}

int *createNonDirichletIndices(int nn,int nd,int *dirich_indices){
    int *ndi = new int[3*nn-nd];
    int pos = 0;
    for(int i=1;i<=3*nn;i++)
        if(!findIndex(i,nd,dirich_indices)){
            ndi[pos] = i;
            pos++;
        }
    return ndi;
}

void writeResults(mesh m,Vector T,char *filename){
    char outputfilename[150];

    int nn = m.getSize(NODES);
    int nd = m.getSize(DIRICHLET);
    int *dirich_indices = m.getDirichletIndices();
    int *non_dirich_indices = createNonDirichletIndices(nn,nd,dirich_indices);
    condition *dirich = m.getDirichlet();
    ofstream file;

    addExtension(outputfilename,filename,".post.res");
    file.open(outputfilename);

    file << "GiD Post Results File 1.0\n";

    file << "Result \"Velocity\" \"Load Case 1\" 1 Vector OnNodes\nComponentNames \"u\" \"v\"\nValues\n";

    for(int i=0;i<nn;i++){
        int d_index = getIndex(i+1,nd,dirich_indices);
        if(d_index != -1)
            file << i+1 << " " << dirich[d_index].getValue() << " ";
        else{
            int T_index = getIndex(i+1,3*nn-nd,non_dirich_indices);
            file << i+1 << " " << T.at(T_index) << " ";
        }
        d_index = getIndex(i+1+nn,nd,dirich_indices);
        if(d_index != -1)
            file << dirich[d_index].getValue() << "\n";
        else{
            int T_index = getIndex(i+1+nn,3*nn-nd,non_dirich_indices);
            file << T.at(T_index) << "\n";
        }
    }

    file << "End values\n";

    file << "\nResult \"Pressure\" \"Load Case 1\" 1 Scalar OnNodes\nComponentNames \"p\"\nValues\n";

    for(int i=0;i<nn;i++){
        int d_index = getIndex(i+1+2*nn,nd,dirich_indices);
        if(d_index != -1)
            file << i+1 << " " << dirich[d_index].getValue() << "\n";
        else{
            int T_index = getIndex(i+1+2*nn,3*nn-nd,non_dirich_indices);
            file << i+1 << " " << T.at(T_index) << "\n";
        }
    }

    file << "End values\n";

    file.close();
}