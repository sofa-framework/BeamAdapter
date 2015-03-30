/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 beta 4      *
*                (c) 2006-2009 MGH, INRIA, USTL, UJF, CNRS                    *
*                                                                             *
* This library is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This library is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this library; if not, write to the Free Software Foundation,     *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.          *
*******************************************************************************
*                               SOFA :: Modules                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
//
//
// Description:
//
//
// Author: Christian Duriez, INRIA
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef SOFA_COMPONENT_ENGINE_WIRERESTSHAPE_INL
#define SOFA_COMPONENT_ENGINE_WIRERESTSHAPE_INL
#define VERIF 1

#include "WireRestShape.h"
#include <sofa/core/behavior/MechanicalState.h>
#include <SofaBaseTopology/EdgeSetGeometryAlgorithms.h>
#include <SofaBaseTopology/QuadSetTopologyModifier.h>

#include <sofa/simulation/common/Node.h>
#include <sofa/simulation/common/TopologyChangeVisitor.h>
#include <sofa/core/topology/Topology.h>

#include <iostream>
#include <fstream>

#define PI 3.14159265
#define EPSILON 0.0000000001

namespace sofa
{

namespace component
{

namespace engine
{

template <class DataTypes>
void WireRestShape<DataTypes>::RotateFrameForAlignX(const Quat &input, Vec3 &x, Quat &output)
{
    x.normalize();
    Vec3 x0=input.inverseRotate(x);

    Real cTheta=x0[0];
     Real theta;
    if (cTheta>(1-EPSILON))
    {
        output = input;
    }
    else
    {
        theta=acos(cTheta);
        // axis of rotation
        Vec3 dw(0,-x0[2],x0[1]);
        dw.normalize();

        // computation of the rotation
        Quat inputRoutput;
        inputRoutput.axisToQuat(dw, theta);

        output=input*inputRoutput;
    }
}

template<class DataTypes>
void WireRestShape<DataTypes>::init()
{

    if (f_printLog.getValue())
        sout<<"WireRestShape begin init"<<sendl;

    if(!procedural.getValue())
    {
        //get the mesh loader
        this->getContext()->get(loader);

        if (!loader)
            serr << "Cannot find a mesh loader. Please insert a MeshObjLoader in the same node" << sendl;
        else
        {
            if (f_printLog.getValue()) sout << "Found a mesh with " << loader->edges.getValue().size() << " edges" << sendl;
            InitFromLoader();
            InitRestConfig();
        }
    }

    //////////////////////////////////////////////
    ////////// get and fill local topology ///////
    //////////////////////////////////////////////

     this->getContext()->get(_topology);
    //_topology = this->getContext()->getMeshTopology();


    if(_topology != NULL)
    {
        if (f_printLog.getValue()) sout<<"found topology named "<< _topology->getName()<<sendl;
    }
    else
        serr << "cannot find topology container" << sendl;

    this->getContext()->get(edgeGeo);
    this->getContext()->get(edgeMod);

    if (edgeGeo == NULL)
        serr << "EdgeSetController has no binding EdgeSetGeometryAlgorithms." << sendl;

    if (edgeMod == NULL)
    {
        serr << "EdgeSetController has no binding EdgeSetTopologyModifier." << sendl;
        edgeSetInNode=false;
    }

    // fill topology :
    _topology->clear();
    _topology->cleanup();
    Real dx = this->length.getValue() / numEdges.getValue();

    // add points
    for ( int i=0; i<numEdges.getValue()+1; i++)
        _topology->addPoint( i*dx, 0, 0);
    // add segments
    for (int i=0; i<numEdges.getValue(); i++)
        _topology->addEdge(i,i+1);

    //// get the possible Topological mapping (with tags)
    const sofa::core::objectmodel::TagSet &tags = this->getTags() ;
    for (core::objectmodel::TagSet::const_iterator it=tags.begin();it!=tags.end();++it)
    {
        std::cerr<<"!!!!!!!!!!!! \n ERROR  : NEED TO FIX line 148 in WireRestShape.inl !!!\n!!!!!!!!"<<std::endl;
        dynamic_cast<core::objectmodel::BaseContext *>(this->getContext())->get( edge2QuadMap , *it, sofa::core::objectmodel::BaseContext::SearchRoot );
    }

    if(!edge2QuadMap)
        serr <<"[WARNING] No Edge2QuadTopologicalMapping map found to propagate the topological change to the topological mapping"<<sendl;

    ////////////////////////////////////////////////////////
    ////////// keyPoint list and Density Assignement ///////
    ////////////////////////////////////////////////////////
    sofa::helper::vector<Real> &keyPointList = (*keyPoints.beginEdit());
    if(!keyPointList.size())
    {
        keyPointList.push_back(0.0);
        if(straightLength.getValue()>= 0.001*this->length.getValue() && straightLength.getValue() <=  0.999*length.getValue())
            keyPointList.push_back(straightLength.getValue());
        keyPointList.push_back(length.getValue());
    }

    keyPoints.endEdit();


    if( density.getValue().size() != keyPointList.size()-1)
    {
        sofa::helper::vector<int> &densityList = (*density.beginEdit());

        if(density.getValue().size() > keyPointList.size()-1 )
            densityList.resize(keyPointList.size()-1);
        else
        {
            densityList.clear();

            if(straightLength.getValue()>= 0.001*this->length.getValue() )
            {
                int numNodes = (int) floor(5.0*straightLength.getValue() / length.getValue() );
                densityList.push_back(numNodes);
            }
            if( straightLength.getValue() <=  0.999*length.getValue())
            {
                int numNodes = (int) floor(20.0*(1.0 - straightLength.getValue() / length.getValue()) );
                densityList.push_back(numNodes);
            }
        }
        density.endEdit();
    }

    if(!numEdgesCollis.getValue().size())
    {
        sofa::helper::vector<int> &densityCol =  (*numEdgesCollis.beginEdit());
        densityCol.resize(keyPointList.size()-1);
        for (unsigned int i=0; i<densityCol.size(); i++)
            densityCol[i] = 20;

        numEdgesCollis.endEdit();
    }

    if (f_printLog.getValue())
        sout<<"WireRestShape end init"<<sendl;

    // Prepare beam sections
    double r 					= this->_radius1.getValue();
    double rInner 				= this->_innerRadius1.getValue();
    this->beamSection1._r 		= r;
    this->beamSection1._rInner 	= rInner;
    this->beamSection1._Iz		= M_PI*(r*r*r*r - rInner*rInner*rInner*rInner)/4.0;
    this->beamSection1._Iy 		= this->beamSection1._Iz ;
    this->beamSection1._J 		= this->beamSection1._Iz + this->beamSection1._Iy;
    this->beamSection1._A 		= M_PI*(r*r - rInner*rInner);
    this->beamSection1._Asy 	= 0.0;
    this->beamSection1._Asz 	= 0.0;

    r 							= this->_radius2.getValue();
    rInner 						= this->_innerRadius2.getValue();
    this->beamSection2._r 		= r;
    this->beamSection2._rInner 	= rInner;
    this->beamSection2._Iz 		= M_PI*(r*r*r*r - rInner*rInner*rInner*rInner)/4.0;
    this->beamSection2._Iy 		= this->beamSection2._Iz ;
    this->beamSection2._J 		= this->beamSection2._Iz + this->beamSection2._Iy;
    this->beamSection2._A 		= M_PI*(r*r - rInner*rInner);
    this->beamSection2._Asy 	= 0.0;
    this->beamSection2._Asz 	= 0.0;
}


template <class DataTypes>
void WireRestShape<DataTypes>::bwdInit()
{
/*
    sofa::core::objectmodel::BaseContext* context = this->getContext();
    core::behavior::MechanicalState<DataTypes> *mState;

    mState = dynamic_cast< core::behavior::MechanicalState<DataTypes> *> (context->getMechanicalState());
    if (!mState)
        serr << "MechanicalStateController has no binding MechanicalState" << sendl;

    VecCoord& x = (*mState->getX());

    Real step=this->length.getValue()/(x.size()-1);

    x[0].clear();
    for (unsigned int i=1; i<x.size(); i++)
    {
        Real x1= step*i;
        Transform global_H_local1;

        this->getRestTransformOnX(global_H_local1, x1);


        x[i].getCenter() = global_H_local1.getOrigin();
        x[i].getOrientation() = global_H_local1.getOrientation();
     }

  */
}


template <class DataTypes>
void WireRestShape<DataTypes>::releaseWirePart(){

    brokenIn2.setValue(true);

    if ( edgeMod == NULL )
    {
        serr<<" no edgeSetModifier in the node -> cannot do the topological change"<<sendl;
        return;
    }
    ///////// remove the edge that is cut //////
    for ( int i=0; i<_topology->getNbPoints(); i++)
    {
        if( _topology->getPX(i) > this->getReleaseCurvAbs() + EPSILON )
        {
            sofa::helper::vector<unsigned int> edge_remove;
            edge_remove.push_back( i-1 );

            std::cout<<"releaseWirePart()  -> remove edge number "<<i<<std::endl;

           edgeMod->removeEdges(edge_remove,false, false); // remove the single edge and do not remove any point...
          // edgeMod->removeEdgesWarning(edge_remove ) ;
          //


           std::cout<<"WireRestShape _topology name="<<_topology->getName()<<" - numEdges ="<<_topology->getNbEdges()<<std::endl;


           // propagate the topological change to the topological mapping //
           if(edge2QuadMap!=NULL)
           {
               //edge2QuadMap->init();
                edge2QuadMap->updateTopologicalMappingTopDown();
                sofa::component::topology::QuadSetTopologyModifier *quadMod;

                edge2QuadMap->getContext()->get(quadMod);


                quadMod->propagateTopologicalChanges();
                //sofa::simulation::Node *node = dynamic_cast<sofa::simulation::Node*> (edge2QuadMap->getContext());
                //sofa::simulation::TopologyChangeVisitor v()

            }


           _topology->resetTopologyChangeList();


           /*
           while ( )
           sofa::simulation::Node *node = static_cast<sofa::simulation::Node*> (this->getContext());
           obj->updateTopologicalMappingTopDown(); // update the specific TopologicalMapping
           */


           return;
        }

    }

    std::cout<<" Wire Part is brokenIn2... should implement a topo change !"<<std::endl;
}


template <class DataTypes>
void WireRestShape<DataTypes>::getSamplingParameters(helper::vector<Real>& xP_noticeable, helper::vector<int>& nbP_density)
{

    xP_noticeable.clear();
    nbP_density.clear();

    if (brokenIn2.getValue())
    {
        for (unsigned int i=0; i<keyPoints.getValue().size(); i++)
        {
            Real x=keyPoints.getValue()[i];
            if( x + EPSILON > this->getReleaseCurvAbs() )
                break;
            xP_noticeable.push_back(x);
            nbP_density.push_back(density.getValue()[i]);
        }
        xP_noticeable.push_back( this->getReleaseCurvAbs());

        std::cout<<"getSamplingParameters brokenIn2 detected - return  xP_noticeable ="<<xP_noticeable<<" and nbP_density ="<<nbP_density<<std::endl;



    }
    else
    {
        xP_noticeable = keyPoints.getValue();
        nbP_density = density.getValue();
    }


/*
    xP_noticeable.push_back(0.0);

    if(brokenIn2.getValue())
    {
        xP_noticeable.push_back(straightLength.getValue());
        int numNodes = (int) floor(5.0*straightLength.getValue() / length.getValue() );
        nbP_density.push_back(numNodes);
        return;
    }

#ifdef VERIF
    // verif:
    if (straightLength.getValue() > this->length.getValue())
    {
        serr<<"straightLength is not <= length => setting straightLength = length"<<sendl;
        straightLength.setValue( length.getValue() );
    }
#endif

    /////// TODO: Alternative Sampling Strategy (i.e depending on the curve of the wirerRestShape)

    if (straightLength.getValue() >= 0.001*this->length.getValue())
    {
        xP_noticeable.push_back(straightLength.getValue());

        int numNodes = (int) floor(5.0*straightLength.getValue() / length.getValue() );
        nbP_density.push_back(numNodes);
    }

    if( straightLength.getValue() <= 0.999*this->length.getValue())
    {
        xP_noticeable.push_back( this->length.getValue() );

        int numNodes = (int) floor(20.0*(1.0 - straightLength.getValue() / length.getValue()) );
        nbP_density.push_back(numNodes);

    }
    */

}

template <class DataTypes>
void WireRestShape<DataTypes>::getCollisionSampling(Real &dx, const Real &x_curv)
{
    unsigned int numLines;
     Real x_used = x_curv - EPSILON;
     if(x_used>length.getValue())
         x_used=length.getValue();

     if(x_used<0.0)
         x_used=0.0;


     // verify that size of numEdgesCollis  =  size of keyPoints-1
     if( numEdgesCollis.getValue().size() != keyPoints.getValue().size()-1)
     {
         serr<<"ooooo\n ooo Problem size of numEdgesCollis ()" << numEdgesCollis.getValue().size() << " !=  size of keyPoints-1 " << keyPoints.getValue().size()-1 <<sendl;
         numLines = (unsigned int)numEdgesCollis.getValue()[0];
         dx=length.getValue()/numLines;
         return;
     }


     for (unsigned int i=1; i<this->keyPoints.getValue().size(); i++)
     {
         if( x_used < this->keyPoints.getValue()[i] )
         {
             numLines = (unsigned int)numEdgesCollis.getValue()[i-1];
             dx=(this->keyPoints.getValue()[i] - this->keyPoints.getValue()[i-1])/numLines;
             return;
         }
     }

     dx=length.getValue()/20;
     serr<<" problem is  getCollisionSampling : x_curv "<<x_used<<" is not between keyPoints"<<keyPoints.getValue()<<std::endl;


}


template <class DataTypes>
void WireRestShape<DataTypes>::getRestTransformOnX(Transform &global_H_local, const Real &x)
{
    Real x_used = x - EPSILON;

    if(x_used>length.getValue())
        x_used=length.getValue();

    if(x_used<0.0)
        x_used=0.0;

    if( x_used < straightLength.getValue())
    {
        global_H_local.set(Vec3(x_used, 0.0, 0.0 ), sofa::defaulttype::Quat());
        return;
    }

    if(procedural.getValue())
    {
        Real projetedLength = spireDiameter.getValue()*PI;
        Real lengthSpire=sqrt(spireHeight.getValue()*spireHeight.getValue() + projetedLength*projetedLength );
            // angle in the z direction
        Real phi= atan(spireHeight.getValue()/projetedLength);

        Quat Qphi;
        Qphi.axisToQuat(Vec3(0,0,1),phi);

        // spire angle (if theta=2*PI, there is a complete spire between startx and x_used)
        Real lengthCurve= x_used-straightLength.getValue();
        Real numSpire=lengthCurve/lengthSpire;
        Real theta= 2*PI*numSpire;
        //std::cout<<"numSpire = "<<numSpire<<"  - theta = "<<theta<<" lengthSpire ="<<lengthSpire<<std::endl;

        // computation of the Quat
        Quat Qtheta;
        Qtheta.axisToQuat(Vec3(0,1,0),theta);
        Quat newSpireQuat = Qtheta*Qphi;


        // computation of the position
        Real radius=spireDiameter.getValue()/2.0;
        //std::cout<<"radius ="<<radius<<std::endl;
        Vec3 PosEndCurve(radius*sin(theta), numSpire*spireHeight.getValue(), radius*(cos(theta)-1)  );

        Vec3 SpirePos=PosEndCurve + Vec3(straightLength.getValue(),0,0);


        global_H_local.set(SpirePos,newSpireQuat);
    }
    else
    {
        x_used = x_used - straightLength.getValue();
        x_used = x_used/(length.getValue()-straightLength.getValue()) * absOfGeometry;

        //std::cout<<"-------------------"<<std::endl;

        Coord p;
        this->getRestPosNonProcedural(x_used,p);
        Vec3 PosEndCurve = p.getCenter();

        //std::cout<<"PosEndCurve"<<PosEndCurve<<std::endl;

        Quat ExtremityQuat = p.getOrientation();

        Vec3 ExtremityPos = PosEndCurve + Vec3(straightLength.getValue(),0,0);

        global_H_local.set(ExtremityPos,ExtremityQuat);
    }
}


template <class DataTypes>
void WireRestShape<DataTypes>::getYoungModulusAtX(const Real& x_curv, Real& youngModulus, Real& cPoisson)
{
    //Initialization
    Real _E1, _E2;
    youngModulus = 0.0;
    cPoisson = 0.0;

    //Get the two possible values of the Young modulus
    _E1 = this->_youngModulus1.getValue();
    _E2 = this->_youngModulus2.getValue();

    //Get User data
    cPoisson = this->_poissonRatio.getValue();

    //Depending on the position of the beam, determine the Young modulus
    if(x_curv <= this->straightLength.getValue())
    {
        youngModulus = _E1;
    }
    else
    {
        if(_E2 == 0.0)
        {
            youngModulus = _E1;
        //	std::cout<<"WARNING : second Young Modulus defined as zero -- only E1 is used"<<std::endl;	// Uncomment if you want a message flood
        }
        else
            youngModulus = _E2;
    }

    return;
}


template <class DataTypes>
void WireRestShape<DataTypes>::getInterpolationParam(const Real& x_curv, Real &_rho, Real &_A, Real &_Iy , Real &_Iz, Real &_Asy, Real &_Asz, Real &_J)
{
    if(x_curv <= this->straightLength.getValue())
    {
        if(_massDensity1.isSet())
            _rho = _massDensity1.getValue();

        if(_radius1.isSet())
        {
            _A		=beamSection1._A;
            _Iy		=beamSection1._Iy;
            _Iz		=beamSection1._Iz;
            _Asy	=beamSection1._Asy;
            _Asz	=beamSection1._Asz;
            _J		=beamSection1._J;
        }
    }
    else
    {
        if(_massDensity2.isSet())
            _rho = _massDensity2.getValue();
        else if(_massDensity1.isSet())
            _rho = _massDensity1.getValue();

        if(_radius2.isSet())
        {
            _A		=beamSection2._A;
            _Iy		=beamSection2._Iy;
            _Iz		=beamSection2._Iz;
            _Asy	=beamSection2._Asy;
            _Asz	=beamSection2._Asz;
            _J		=beamSection2._J;
        }
        else if(_radius1.isSet())
        {
            _A		=beamSection1._A;
            _Iy		=beamSection1._Iy;
            _Iz		=beamSection1._Iz;
            _Asy	=beamSection1._Asy;
            _Asz	=beamSection1._Asz;
            _J		=beamSection1._J;
        }
    }
}

template <class DataTypes>
bool WireRestShape<DataTypes>::checkTopology()
{
    if (!loader->edges.getValue().size())
    {
        serr << "There is no edges in the topology loaded by " << loader->getName() << sendl;
        return false;
    }

    if (loader->triangles.getValue().size())
    {
        serr << "There are triangles in the topology loaded by " << loader->getName() << sendl;
        return false;
    }

    if (loader->quads.getValue().size())
    {
        serr << "There are quads in the topology loaded by " << loader->getName() << sendl;
        return false;
    }

    if (loader->polygons.getValue().size())
    {
        serr << "There are polygons in the topology loaded by " << loader->getName() << sendl;
        return false;
    }

    /// \todo check if the topology is like a wire


    return true;
}



template <class DataTypes>
void WireRestShape<DataTypes>::InitFromLoader()
{
    if (!checkTopology())
        return;

    sofa::helper::vector<Vec3> vertices;
    sofa::helper::vector<Vec2> edges;

    //get the topology position
    typedef  sofa::helper::vector<sofa::defaulttype::Vec<3,SReal> > topoPosition;
    topoPosition &topoVertices = (*loader->positions.beginEdit());

    //copy the topology edges in a local vector
    typedef  sofa::helper::vector<sofa::core::topology::Topology::Edge > topoEdge;
    topoEdge &topoEdges = (*loader->edges.beginEdit());
    for (topoEdge::iterator it = topoEdges.begin(); it < topoEdges.end(); it++)
        edges.push_back(Vec2((*it)[0], (*it)[1]));
    loader->edges.endEdit();

    /** renumber the vertices  **/
   sofa::helper::vector<unsigned int> verticesConnexion; //gives the number of edges connected to a vertex
   for(unsigned int i =0; i < topoVertices.size(); i++)
       verticesConnexion.push_back(2);

   for(unsigned int i = 0; i < edges.size(); i++)
   {
        Vec2 ed = edges[i];
        unsigned int e1 = floor(ed[0]);
        unsigned int e2 = floor(ed[1]);
        verticesConnexion[e1]--;
        verticesConnexion[e2]--;
   }
   if (this->f_printLog.getValue())
       sout << "Successfully compute the vertex connexion" << sendl;

   // check for the first corner of the edge
   unsigned int firstIndex = 0;
   bool found = false;
   while((firstIndex < verticesConnexion.size()) && !found)
   {
       if(verticesConnexion[firstIndex] == 1)
           found = true;
       else
           firstIndex++;
   }

   if(firstIndex == verticesConnexion.size())
   {
       serr << "The first vertex of the beam structure is not found, probably because of a closed structure" << sendl;
       return;
   }

   vertices.push_back(topoVertices[firstIndex]);

    while(edges.size() > 0)
    {
        vecIt it = edges.begin();
        vecIt end = edges.end();

        bool notFound = true;
        while (notFound && (it != end))
        {
            Vec2 ed = (*it);
            vecIt toDel = it;
            it++;
            if(ed[0] == firstIndex)
            {
                vertices.push_back(topoVertices[ed[1]]);
                firstIndex = ed[1];
                //std::cout << firstIndex << " added " << std::endl;
                edges.erase(toDel);
                notFound = false;

            }
            else if(ed[1] == firstIndex)
            {
                vertices.push_back(topoVertices[ed[0]]);
                firstIndex = ed[0];
                //std::cout << firstIndex << " added " << std::endl;
                edges.erase(toDel);
                notFound = false;
            }
        }
    }

    if (this->f_printLog.getValue())
        sout << "Successfully computed the topology" << sendl;

    localRestPositions = vertices;

    for(unsigned int i = 0; i < localRestPositions.size() - 1; i++)
        localRestPositions[i] *= NonProceduralScale.getValue();

    loader->positions.endEdit();
}


template <class DataTypes>
void WireRestShape<DataTypes>::InitRestConfig()
{
    curvAbs.clear();
    double tot = 0;
    curvAbs.push_back(0);
    Quat input, output;
    input.identity();
    localRestTransforms.resize(localRestPositions.size());
    localRestTransforms[0].setOrigin(Vec3(0,0,0));
    localRestTransforms[0].setOrientation(input);

    for(unsigned int i = 0; i < localRestPositions.size() - 1; i++)
    {
        Vec3 vec = localRestPositions[i+1] - localRestPositions[i];
        double norm = vec.norm();
        tot += norm;

        this->RotateFrameForAlignX(input, vec, output);

        input = output;

        localRestTransforms[i+1].setOrientation(output);

        Vec3 localPos = localRestPositions[i+1] - localRestPositions[0];

        localRestTransforms[i+1].setOrigin(localPos);

        //        if (this->f_printLog.getValue())
        //        	sout <<"localRestTransforms ="<<localRestTransforms[i]<<sendl;

        curvAbs.push_back(tot);

        //std::cout << "++++++++++++----------------++++++++++++++++++localRestPositions = " << localRestPositions[i] << std::endl;
        //std::cout << "++++++++++++----------------++++++++++++++++++curveAbs = " << tot << std::endl;
    }
    absOfGeometry = tot;

    Real newLength = straightLength.getValue() + absOfGeometry;
    length.setValue(newLength);

    if (f_printLog.getValue())
        sout <<"Length of the loaded shape = "<< absOfGeometry << ", total length with straight length = " << newLength << sendl;
}


template <class DataTypes>
void WireRestShape<DataTypes>::getRestPosNonProcedural(Real& abs, Coord &p)
{
   /*** find the range which includes the "requested" abs ***/
   double startingAbs = 0; unsigned int index = 0;

   while ((startingAbs < abs) && (index < localRestPositions.size()))
   {
       index++;
       startingAbs = curvAbs[index];
    }

   //std::cout<<"index/localRestPositionsMAX = "<<index<<"/"<<localRestPositions.size()<<std::endl;
   //std::cout<<"abs = "<<abs<<std::endl;

   /*** OOB ***/
   if(abs > startingAbs)
   {
       serr << "abs = "<<abs<<" et startingAbs = "<< startingAbs<< sendl;
       serr << "[Warning] Out of bound position request" << sendl;
       return ;
   }
   else /*** Expected case ***/
   {
       /*** restAbs[index-1]   abs                       restAbs[index]
                    |-----------|----------------------------|
                       alpha         one_minu_alpha
       ***/


       Real alpha, one_minus_alpha;
       Vec3 result;

      /// std::cout<<" abs ="<<abs<<" - curvAbs["<<index-1<<"]:"<<curvAbs[index-1]<<" - curvAbs["<<index<<"]:"<<curvAbs[index]<<" - curvAbs["<<index+1<<"]:"<<curvAbs[index+1]<<std::endl;

       alpha = (abs - curvAbs[index-1] ) / (curvAbs[index] - curvAbs[index-1]);
       one_minus_alpha = 1 - alpha;
       result = localRestTransforms[index - 1].getOrigin() * one_minus_alpha + localRestTransforms[index].getOrigin() * alpha;
       Quat slerp;
       slerp.slerp( localRestTransforms[index - 1].getOrientation(),  localRestTransforms[index].getOrientation(), alpha, true );

       slerp.normalize();

       p.getCenter() = result;

       p.getOrientation() = slerp;
   }
}



template <class DataTypes>
void WireRestShape<DataTypes>::computeOrientation(const Vec3& AB, const Quat& Q, Quat &result)
{
    Vec3 PQ = AB;
    Quat quat = Q;

    Vec3 x = quat.rotate(Vec3(1,0,0));
    PQ.normalize();

    if (dot(x, PQ) > 0.9999999)
        result = Q;

    Vec3 y;
    double alpha;

    if (dot(x, PQ) < -0.9999999)
    {
        y = quat.rotate(Vec3(0,0,1));
        alpha = M_PI;
    }
    else
    {
        y = cross(x, PQ);
        y.normalize();
        alpha = acos(dot(x, PQ));
    }

    Quat qaux = Quat(y, alpha);
    result = qaux * quat;

}

template<class DataTypes>
void WireRestShape<DataTypes>::draw(const core::visual::VisualParams* /*vparams*/)
{
    if (!drawRestShape.getValue())
        return;

    glDisable(GL_LIGHTING);
    glColor3d(1.0,0.0,0.0);
    glPointSize(10.0);
    for (unsigned int i = 0 ; i < localRestPositions.size(); i++)
    {
        glBegin(GL_POINTS);
            glVertex3d(localRestPositions[i][0],localRestPositions[i][1],localRestPositions[i][2]);
        glEnd();
    }
    glPointSize(1.0);
    glEnable(GL_LIGHTING);

}



}// namespace engine


} // namespace component

} // namespace sofa

#endif /* SOFA_COMPONENT_ENGINE_WIRERESTSHAPE_INL */
