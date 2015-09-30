#ifndef vertex_H
#define vertex_H

#include "helper.h"

using namespace std;


int getLeg(int Vertexid, int leg);

int flipLeg(int Vertexid, int leg);

void getAllLegs(int Vertexid, int& l0, int& l1, int& l2, int& l3);

float getDiagonalOffset(const float& J12, const float& h1, const float& h2);

float getDiagonalWeight(const float& J12, const float& h1, const float& h2, int Vertexid);

int getCompatibleDiagVertex(int s0, int s1);

int getVertexType(int Vertexid);

int isFlippedLeg(int Vertexid, int leg);

int getFlippedLeg(int Vertexid);

void getVertices(const float& _J12, const float& _h1, const float& _h2, const float& _d1, const float _d2 , float& diagshift, array<float,16>& VWeights, array<float,4>& diagVWeights);

int getVertexID(int s0, int s1, int s2, int s3);

int getVertexID(int s0, int s1, int otype);

int VertexToOperator(int vertexID);

int FlipVertex(int vID, int enleg, int exleg);


#endif
