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

#include <sofa/core/behavior/MechanicalState.h>
#include <SofaBaseTopology/QuadSetTopologyModifier.h>
#include <SofaBaseTopology/EdgeSetGeometryAlgorithms.h>
#include <SofaTopologyMapping/Edge2QuadTopologicalMapping.h>
#include <SofaLoader/MeshObjLoader.h>

#include <sofa/simulation/Node.h>
#include <sofa/simulation/TopologyChangeVisitor.h>
#include <sofa/core/topology/Topology.h>
#include <sofa/core/visual/VisualParams.h>

#include <sofa/helper/gl/template.h>
//#include <sofa/helper/system/glut.h>

#include "WireRestShape.h"

#include <iostream>
#include <fstream>

#define EPSILON 0.0000000001
#define VERIF 1

namespace sofa
{

namespace component
{

namespace engine
{

namespace _wirerestshape_
{

using sofa::helper::vector ;
using sofa::core::objectmodel::TagSet ;
using sofa::core::objectmodel::BaseContext ;

/*!
 * @brief Default Constructor.
 */
template <class DataTypes>
WireRestShape<DataTypes>::WireRestShape() :
    //TODO(dmarchal 2017-05-17) not sure that procedural & nonProceduralScale are very understandable name...are they exclusives ?
    //if so have look in my comment in the init section.
    d_isAProceduralShape( initData(&d_isAProceduralShape,(bool)true,"isAProceduralShape","is the guidewire shape mathemetically defined ?") )
  , d_nonProceduralScale( initData ( &d_nonProceduralScale, (Real)1.0, "nonProceduralScale", "scale of the model defined by file" ) )
  , d_length(initData(&d_length, (Real)1.0, "length", "total length of the wire instrument"))
  , d_straightLength(initData(&d_straightLength, (Real)0.0, "straightLength", "length of the initial straight shape"))
  , d_spireDiameter(initData(&d_spireDiameter, (Real)0.1, "spireDiameter", "diameter of the spire"))
  , d_spireHeight(initData(&d_spireHeight, (Real)0.01, "spireHeight", "height between each spire"))
  , d_density(initData(&d_density, "densityOfBeams", "density of beams between key points"))
  , d_keyPoints(initData(&d_keyPoints,"keyPoints","key points of the shape (curv absc)"))
  , d_numEdges(initData(&d_numEdges, 10, "numEdges","number of Edges for the visual model"))
  , d_numEdgesCollis(initData(&d_numEdgesCollis,"numEdgesCollis", "number of Edges for the collision model" ))
  , d_poissonRatio(initData(&d_poissonRatio,(Real)0.49,"poissonRatio","Poisson Ratio"))
  , d_youngModulus1(initData(&d_youngModulus1,(Real)5000,"youngModulus","Young Modulus"))
  , d_youngModulus2(initData(&d_youngModulus2,(Real)3000,"youngModulusExtremity","youngModulus for beams at the extremity\nonly if not straight"))
  , d_radius1(initData(&d_radius1,(Real)1.0f,"radius","radius"))
  , d_radius2(initData(&d_radius2,(Real)1.0f,"radiusExtremity","radius for beams at the extremity\nonly if not straight"))
  , d_innerRadius1(initData(&d_innerRadius1,(Real)0.0f,"innerRadius","inner radius if it applies"))
  , d_innerRadius2(initData(&d_innerRadius2,(Real)0.0f,"innerRadiusExtremity","inner radius for beams at the extremity\nonly if not straight"))
  , d_massDensity1(initData(&d_massDensity1,(Real)1.0,"massDensity", "Density of the mass (usually in kg/m^3)" ))
  , d_massDensity2(initData(&d_massDensity2,(Real)1.0,"massDensityExtremity", "Density of the mass at the extremity\nonly if not straight" ))
  , d_brokenIn2(initData(&d_brokenIn2, (bool)false, "brokenIn2", ""))
  , d_drawRestShape(initData(&d_drawRestShape, (bool)false, "draw", "draw rest shape"))
  , edge2QuadMap(nullptr)
{
    d_spireDiameter.setGroup("Procedural");
    d_spireHeight.setGroup("Procedural");
}

template <class DataTypes>
void WireRestShape<DataTypes>::rotateFrameForAlignX(const Quat &input, Vec3 &x, Quat &output)
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
void WireRestShape<DataTypes>::parse(BaseObjectDescription* args)
{
    const char* arg = args->getAttribute("procedural") ;
    if(arg)
    {
        msg_info() << "The attribute 'procedural' has been renamed into 'isAProceduralShape'. " << msgendl
                   << "To remove this warning you need to update your scene and replace 'procedural' with 'isAProceduralShape'" ;

        args->removeAttribute("procedural") ;
        args->setAttribute("isAProceduralShape", arg) ;
    }

    Inherit1::parse(args) ;
}

template<class DataTypes>
void WireRestShape<DataTypes>::init()
{
    /// Have to add because was remove from the .h due to forward déclarations
    edgeGeo = nullptr;
    edge2QuadMap = nullptr;
    loader = nullptr;

    if(!d_isAProceduralShape.getValue())
    {
        /// get the mesh loader
        this->getContext()->get(loader);

        if (!loader)
            msg_error() << "Cannot find a mesh loader. Please insert a MeshObjLoader in the same node" ;
        else
        {
            msg_info() << "Found a mesh with " << loader->d_edges.getValue().size() << " edges" ;
            initFromLoader();
            initRestConfig();
        }
    }

    //////////////////////////////////////////////
    ////////// get and fill local topology ///////
    //////////////////////////////////////////////
    this->getContext()->get(_topology);

    if(_topology != NULL){
        msg_info() << "found topology named "<< _topology->getName() ;
    }else{
        //TODO(dmarchal): Explain here what will happen to the user and how he can remove this message"
        msg_error() << "cannot find topology container";
    }

    this->getContext()->get(edgeGeo);
    this->getContext()->get(edgeMod);

    //TODO(dmarchal): Explain here what will happen to the user and how he can remove this message"
    msg_error_when(edgeGeo == NULL)
            << "EdgeSetController has no binding EdgeSetGeometryAlgorithms." ;

    if (edgeMod == NULL)
    {
        //TODO(dmarchal): Explain here what will happen to the user and how he can remove this message"
        msg_error() << "EdgeSetController has no binding EdgeSetTopologyModifier." ;
        edgeSetInNode=false;
    }

    /// fill topology :
    _topology->clear();
    _topology->cleanup();
    Real dx = this->d_length.getValue() / d_numEdges.getValue();

    /// add points
    for ( int i=0; i<d_numEdges.getValue()+1; i++)
        _topology->addPoint( i*dx, 0, 0);

    /// add segments
    for (int i=0; i<d_numEdges.getValue(); i++)
        _topology->addEdge(i,i+1);

    //// get the possible Topological mapping (with tags)
    const TagSet &tags = this->getTags() ;

    //TODO(dmarchal): replace with for each loop.
    for (TagSet::const_iterator it=tags.begin();it!=tags.end();++it)
    {
        dmsg_error() <<" NEED TO FIX line 148 in WireRestShape.inl " ;
        dynamic_cast<BaseContext *>(this->getContext())->get( edge2QuadMap , *it, BaseContext::SearchRoot );
    }


    //TODO(dmarchal): Explain here what will happen to the user and how he can remove this message"
    if(!edge2QuadMap)
        msg_error()<< "No Edge2QuadTopologicalMapping map found to propagate the topological change to the topological mapping"<<sendl;


    ////////////////////////////////////////////////////////
    ////////// keyPoint list and Density Assignement ///////
    ////////////////////////////////////////////////////////
    vector<Real> &keyPointList = (*d_keyPoints.beginEdit());
    if(!keyPointList.size())
    {
        keyPointList.push_back(0.0);
        if(d_straightLength.getValue()>= 0.001*this->d_length.getValue() && d_straightLength.getValue() <=  0.999*d_length.getValue())
            keyPointList.push_back(d_straightLength.getValue());
        keyPointList.push_back(d_length.getValue());
    }

    d_keyPoints.endEdit();


    if( d_density.getValue().size() != keyPointList.size()-1)
    {
        vector<int> &densityList = (*d_density.beginEdit());

        if(d_density.getValue().size() > keyPointList.size()-1 )
            densityList.resize(keyPointList.size()-1);
        else
        {
            densityList.clear();

            if(d_straightLength.getValue()>= 0.001*this->d_length.getValue() )
            {
                int numNodes = (int) floor(5.0*d_straightLength.getValue() / d_length.getValue() );
                densityList.push_back(numNodes);
            }
            if( d_straightLength.getValue() <=  0.999*d_length.getValue())
            {
                int numNodes = (int) floor(20.0*(1.0 - d_straightLength.getValue() / d_length.getValue()) );
                densityList.push_back(numNodes);
            }
        }
        d_density.endEdit();
    }

    if(!d_numEdgesCollis.getValue().size())
    {
        vector<int> &densityCol =  (*d_numEdgesCollis.beginEdit());
        densityCol.resize(keyPointList.size()-1);
        for (unsigned int i=0; i<densityCol.size(); i++)
            densityCol[i] = 20;

        d_numEdgesCollis.endEdit();
    }

    msg_info() <<"WireRestShape end init" ;

    // Prepare beam sections
    double r 					= this->d_radius1.getValue();
    double rInner 				= this->d_innerRadius1.getValue();
    this->beamSection1._r 		= r;
    this->beamSection1._rInner 	= rInner;
    this->beamSection1._Iz		= M_PI*(r*r*r*r - rInner*rInner*rInner*rInner)/4.0;
    this->beamSection1._Iy 		= this->beamSection1._Iz ;
    this->beamSection1._J 		= this->beamSection1._Iz + this->beamSection1._Iy;
    this->beamSection1._A 		= M_PI*(r*r - rInner*rInner);
    this->beamSection1._Asy 	= 0.0;
    this->beamSection1._Asz 	= 0.0;

    r 							= this->d_radius2.getValue();
    rInner 						= this->d_innerRadius2.getValue();
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
}

template <class DataTypes>
void WireRestShape<DataTypes>::releaseWirePart(){

    d_brokenIn2.setValue(true);

    if ( edgeMod == NULL )
    {
        msg_error() << "no edgeSetModifier in the node -> cannot do the topological change"<<sendl;
        return;
    }
    ///////// remove the edge that is cut //////
    for ( int i=0; i<_topology->getNbPoints(); i++)
    {
        if( _topology->getPX(i) > this->getReleaseCurvAbs() + EPSILON )
        {
            vector<unsigned int> edge_remove;
            edge_remove.push_back( i-1 );

            msg_info() << "releaseWirePart()  -> remove edge number "<< i ;

            edgeMod->removeEdges(edge_remove,false, false); // remove the single edge and do not remove any point...

            msg_info() << "WireRestShape _topology name="<<_topology->getName()<<" - numEdges ="<<_topology->getNbEdges() ;

            // propagate the topological change to the topological mapping //
            if(edge2QuadMap!=NULL)
            {
                edge2QuadMap->updateTopologicalMappingTopDown();
                sofa::component::topology::QuadSetTopologyModifier *quadMod;
                edge2QuadMap->getContext()->get(quadMod);
                quadMod->propagateTopologicalChanges();
            }


            _topology->resetTopologyChangeList();

            return;
        }
    }

    dmsg_info() <<" Wire Part is brokenIn2... should implement a topo change !" ;
}


template <class DataTypes>
void WireRestShape<DataTypes>::getSamplingParameters(vector<Real>& xP_noticeable,
                                                     vector<int>& nbP_density) const
{

    xP_noticeable.clear();
    nbP_density.clear();

    if (d_brokenIn2.getValue())
    {
        for (unsigned int i=0; i<d_keyPoints.getValue().size(); i++)
        {
            Real x=d_keyPoints.getValue()[i];
            if( x + EPSILON > getReleaseCurvAbs() )
                break;
            xP_noticeable.push_back(x);
            nbP_density.push_back(d_density.getValue()[i]);
        }
        xP_noticeable.push_back( getReleaseCurvAbs());

        dmsg_info() <<"getSamplingParameters brokenIn2 detected - return  xP_noticeable ="<<xP_noticeable<<" and nbP_density ="<<nbP_density ;
    }
    else
    {
        xP_noticeable = d_keyPoints.getValue();
        nbP_density = d_density.getValue();
    }
}

template <class DataTypes>
void WireRestShape<DataTypes>::getCollisionSampling(Real &dx, const Real &x_curv)
{
    unsigned int numLines;
    Real x_used = x_curv - EPSILON;
    if(x_used>d_length.getValue())
        x_used=d_length.getValue();

    if(x_used<0.0)
        x_used=0.0;

    // verify that size of numEdgesCollis  =  size of keyPoints-1
    if( d_numEdgesCollis.getValue().size() != d_keyPoints.getValue().size()-1)
    {
        msg_error() << "Problem size of numEdgesCollis ()" << d_numEdgesCollis.getValue().size() << " !=  size of keyPoints-1 " << d_keyPoints.getValue().size()-1 ;
        numLines = (unsigned int)d_numEdgesCollis.getValue()[0];
        dx=d_length.getValue()/numLines;
        return;
    }


    for (unsigned int i=1; i<this->d_keyPoints.getValue().size(); i++)
    {
        if( x_used < this->d_keyPoints.getValue()[i] )
        {
            numLines = (unsigned int)d_numEdgesCollis.getValue()[i-1];
            dx=(this->d_keyPoints.getValue()[i] - this->d_keyPoints.getValue()[i-1])/numLines;
            return;
        }
    }

    dx=d_length.getValue()/20;
    msg_error() << " problem is  getCollisionSampling : x_curv "<<x_used<<" is not between keyPoints"<<d_keyPoints.getValue() ;
}


template <class DataTypes>
void WireRestShape<DataTypes>::getRestTransformOnX(Transform &global_H_local, const Real &x)
{
    Real x_used = x - EPSILON;

    if(x_used>d_length.getValue())
        x_used=d_length.getValue();

    if(x_used<0.0)
        x_used=0.0;

    if( x_used < d_straightLength.getValue())
    {
        global_H_local.set(Vec3(x_used, 0.0, 0.0 ), Quat());
        return;
    }

    if(d_isAProceduralShape.getValue())
    {
        Real projetedLength = d_spireDiameter.getValue()*M_PI;
        Real lengthSpire=sqrt(d_spireHeight.getValue()*d_spireHeight.getValue() + projetedLength*projetedLength );
        // angle in the z direction
        Real phi= atan(d_spireHeight.getValue()/projetedLength);

        Quat Qphi;
        Qphi.axisToQuat(Vec3(0,0,1),phi);

        // spire angle (if theta=2*PI, there is a complete spire between startx and x_used)
        Real lengthCurve= x_used-d_straightLength.getValue();
        Real numSpire=lengthCurve/lengthSpire;
        Real theta= 2*M_PI*numSpire;

        // computation of the Quat
        Quat Qtheta;
        Qtheta.axisToQuat(Vec3(0,1,0),theta);
        Quat newSpireQuat = Qtheta*Qphi;


        // computation of the position
        Real radius=d_spireDiameter.getValue()/2.0;
        Vec3 PosEndCurve(radius*sin(theta), numSpire*d_spireHeight.getValue(), radius*(cos(theta)-1)  );
        Vec3 SpirePos=PosEndCurve + Vec3(d_straightLength.getValue(),0,0);

        global_H_local.set(SpirePos,newSpireQuat);
    }
    else
    {
        x_used = x_used - d_straightLength.getValue();
        x_used = x_used/(d_length.getValue()-d_straightLength.getValue()) * m_absOfGeometry;

        Coord p;
        this->getRestPosNonProcedural(x_used,p);
        Vec3 PosEndCurve = p.getCenter();
        Quat ExtremityQuat = p.getOrientation();
        Vec3 ExtremityPos = PosEndCurve + Vec3(d_straightLength.getValue(),0,0);

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
    _E1 = this->d_youngModulus1.getValue();
    _E2 = this->d_youngModulus2.getValue();

    //Get User data
    cPoisson = this->d_poissonRatio.getValue();

    //Depending on the position of the beam, determine the Young modulus
    if(x_curv <= this->d_straightLength.getValue())
    {
        double varyingYoungModulus = _E1*(d_straightLength.getValue()-x_curv)/d_straightLength.getValue() + _E1/20.0;

        youngModulus = _E1;
        msg_warning() << "youngModulus " <<  youngModulus;
        msg_warning() << "varyingYoungModulus " <<  varyingYoungModulus;
        // Electrode array parameters
        double softTipLength = 0.3;
        double lengthElectrode = 0.5;
        double lengthSilicone = 0.7;
        double arrayLength = 25;
        double lengthSingleElectrodeAndSiliconeSection = lengthElectrode+lengthSilicone;
        double lengthFromTip = d_straightLength.getValue()-x_curv;
        msg_warning() << "x_curv is : " << x_curv;
        msg_warning() << "lengthFromTip is : " << lengthFromTip;
        if (lengthFromTip > arrayLength) // longer than array itself so make it quasi-rigid to mimic surgeon pushing
        {
            youngModulus = _E1*1000;
        }
        else if (lengthFromTip > softTipLength) // main body of the array
        {
            double caracLength = (lengthFromTip - softTipLength) - floor((lengthFromTip - softTipLength)/lengthSingleElectrodeAndSiliconeSection)*lengthSingleElectrodeAndSiliconeSection;
            if (caracLength<0.5) // we are at an electrode
            {
                youngModulus = _E1*100;
                msg_warning() << "RIGID: " << x_curv;
                msg_warning() << "youngModulus: " << youngModulus;

            }
            else // we are at a silicone section
            {
                youngModulus = varyingYoungModulus;
                msg_warning() << "SOFT: " << x_curv;
                msg_warning() << "youngModulus: " << youngModulus;

            }
        }     
        else  // this is the very tip of the array
        {
           youngModulus = varyingYoungModulus;
        }
    }
    else
    {
        if(_E2 == 0.0)
            youngModulus = _E1;
        else
            youngModulus = _E2;
    }
    return;
}


template <class DataTypes>
void WireRestShape<DataTypes>::getInterpolationParam(const Real& x_curv, Real &_rho, Real &_A, Real &_Iy , Real &_Iz, Real &_Asy, Real &_Asz, Real &_J)
{
    if(x_curv <= this->d_straightLength.getValue())
    {
        if(d_massDensity1.isSet())
            _rho = d_massDensity1.getValue();

        if(d_radius1.isSet())
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
        if(d_massDensity2.isSet())
            _rho = d_massDensity2.getValue();
        else if(d_massDensity1.isSet())
            _rho = d_massDensity1.getValue();

        if(d_radius2.isSet())
        {
            _A		=beamSection2._A;
            _Iy		=beamSection2._Iy;
            _Iz		=beamSection2._Iz;
            _Asy	=beamSection2._Asy;
            _Asz	=beamSection2._Asz;
            _J		=beamSection2._J;
        }
        else if(d_radius1.isSet())
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
    if (!loader->d_edges.getValue().size())
    {
        msg_error() << "There is no edges in the topology loaded by " << loader->getName() ;
        return false;
    }

    if (loader->d_triangles.getValue().size())
    {
        msg_error() << "There are triangles in the topology loaded by " << loader->getName() ;
        return false;
    }

    if (loader->d_quads.getValue().size())
    {
        msg_error() << "There are quads in the topology loaded by " << loader->getName() ;
        return false;
    }

    if (loader->d_polygons.getValue().size())
    {
        msg_error() << "There are polygons in the topology loaded by " << loader->getName() ;
        return false;
    }

    //TODO(dmarchal 2017-05-17) when writing a TODO please specify:
    // who will do that
    // when it will be done
    /// \todo check if the topology is like a wire


    return true;
}



template <class DataTypes>
void WireRestShape<DataTypes>::initFromLoader()
{
    if (!checkTopology())
        return;

    vector<Vec3> vertices;
    vector<Vec2> edges;

    //get the topology position
    typedef  vector<sofa::defaulttype::Vec<3,SReal> > topoPosition;
    topoPosition &topoVertices = (*loader->d_positions.beginEdit());

    //copy the topology edges in a local vector
    typedef  vector<sofa::core::topology::Topology::Edge > topoEdge;
    topoEdge &topoEdges = (*loader->d_edges.beginEdit());
    for (topoEdge::iterator it = topoEdges.begin(); it < topoEdges.end(); it++)
        edges.push_back(Vec2((*it)[0], (*it)[1]));
    loader->d_edges.endEdit();

    /** renumber the vertices  **/
    vector<unsigned int> verticesConnexion; //gives the number of edges connected to a vertex
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

    msg_info() << "Successfully compute the vertex connexion" ;

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
        msg_error() << "The first vertex of the beam structure is not found, probably because of a closed structure" ;
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
                edges.erase(toDel);
                notFound = false;

            }
            else if(ed[1] == firstIndex)
            {
                vertices.push_back(topoVertices[ed[0]]);
                firstIndex = ed[0];
                edges.erase(toDel);
                notFound = false;
            }
        }
    }

    msg_info() << "Successfully computed the topology" ;

    m_localRestPositions = vertices;

    for(unsigned int i = 0; i < m_localRestPositions.size() - 1; i++)
        m_localRestPositions[i] *= d_nonProceduralScale.getValue();

    loader->d_positions.endEdit();
}


template <class DataTypes>
void WireRestShape<DataTypes>::initRestConfig()
{
    m_curvAbs.clear();
    double tot = 0;
    m_curvAbs.push_back(0);
    Quat input, output;
    input.identity();
    m_localRestTransforms.resize(m_localRestPositions.size());
    m_localRestTransforms[0].setOrigin(Vec3(0,0,0));
    m_localRestTransforms[0].setOrientation(input);

    for(unsigned int i = 0; i < m_localRestPositions.size() - 1; i++)
    {
        Vec3 vec = m_localRestPositions[i+1] - m_localRestPositions[i];
        double norm = vec.norm();
        tot += norm;

        this->rotateFrameForAlignX(input, vec, output);

        input = output;

        m_localRestTransforms[i+1].setOrientation(output);

        Vec3 localPos = m_localRestPositions[i+1] - m_localRestPositions[0];

        m_localRestTransforms[i+1].setOrigin(localPos);

        m_curvAbs.push_back(tot);
    }
    m_absOfGeometry = tot;

    Real newLength = d_straightLength.getValue() + m_absOfGeometry;
    d_length.setValue(newLength);

    msg_info() <<"Length of the loaded shape = "<< m_absOfGeometry << ", total length with straight length = " << newLength ;
}


template <class DataTypes>
void WireRestShape<DataTypes>::getRestPosNonProcedural(Real& abs, Coord &p)
{
    /*** find the range which includes the "requested" abs ***/
    double startingAbs = 0; unsigned int index = 0;

    while ((startingAbs < abs) && (index < m_localRestPositions.size()))
    {
        index++;
        startingAbs = m_curvAbs[index];
    }

    /*** OOB ***/
    if(abs > startingAbs)
    {
        msg_error() << "abs = "<<abs<<" et startingAbs = "<< startingAbs<< msgendl
                    << "Out of bound position request" ;
        return ;
    }
    else /*** Expected case ***/
    {
        Real alpha, one_minus_alpha;
        Vec3 result;

        alpha = (abs - m_curvAbs[index-1] ) / (m_curvAbs[index] - m_curvAbs[index-1]);
        one_minus_alpha = 1 - alpha;
        result = m_localRestTransforms[index - 1].getOrigin() * one_minus_alpha + m_localRestTransforms[index].getOrigin() * alpha;
        Quat slerp;
        slerp.slerp( m_localRestTransforms[index - 1].getOrientation(),  m_localRestTransforms[index].getOrientation(), alpha, true );

        slerp.normalize();

        p.getCenter() = result;

        p.getOrientation() = slerp;
    }
}

template <class DataTypes>
typename WireRestShape<DataTypes>::Real WireRestShape<DataTypes>::getLength()
{
    if(d_brokenIn2.getValue())
        return d_straightLength.getValue();
    else
        return d_length.getValue();
}

template <class DataTypes>
void WireRestShape<DataTypes>::getNumberOfCollisionSegment(Real &dx, unsigned int &numLines)
{
    numLines = 0;
    for (unsigned i=0; i<d_numEdgesCollis.getValue().size(); i++)
    {
        numLines += (unsigned int)d_numEdgesCollis.getValue()[i];
    }
    dx=d_length.getValue()/numLines;
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
    if (!d_drawRestShape.getValue())
        return;

    glDisable(GL_LIGHTING);
    glColor3d(1.0,0.0,0.0);
    glPointSize(10.0);
    for (unsigned int i = 0 ; i < m_localRestPositions.size(); i++)
    {
        glBegin(GL_POINTS);
        glVertex3d(m_localRestPositions[i][0],m_localRestPositions[i][1],m_localRestPositions[i][2]);
        glEnd();
    }
    glPointSize(1.0);
    glEnable(GL_LIGHTING);

}

} // namespace _wirerestshape_
using _wirerestshape_::WireRestShape;

} // namespace engine

} // namespace component

} // namespace sofa

#endif /* SOFA_COMPONENT_ENGINE_WIRERESTSHAPE_INL */
