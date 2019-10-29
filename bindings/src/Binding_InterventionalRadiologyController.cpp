#include <pybind11/pybind11.h>
#include "../../component/controller/InterventionalRadiologyController.h"
#include "Binding_InterventionalRadiologyController.h"
#include "Binding_InterventionalRadiologyController_doc.h"



namespace sofapython3
{
namespace py { using namespace pybind11; }
using sofa::core::objectmodel::BaseObject;
using sofa::component::controller::InterventionalRadiologyController;

bool modifyTopology(InterventionalRadiologyController<Rigid3Types>* self){
    bool topoModified;
    topoModified = self->modifyTopology();
    return topoModified;    
}

// void fixFirstNodesWithUntil(InterventionalRadiologyController<Rigid3Types>* self){
//     
//     unsigned int firstSimulatedNode;
//     
//     self->fixFirstNodesWithUntil(firstSimulatedNode);
// }

void moduleAddInterventionalRadiologyController(py::module &m)
{

    py::class_<InterventionalRadiologyController<Rigid3Types>,
            sofa::core::behavior::BaseController,
            py_shared_ptr<InterventionalRadiologyController<Rigid3Types>>> p(m, "InterventionalRadiology");

    p.def("getTotalNbEdges",&InterventionalRadiologyController<Rigid3Types>::getTotalNbEdges,
          sofapython3::doc::controller::InterventionalRadiologyController::InterventionalRadiologyControllerClass);


    p.def("modifyTopology", modifyTopology);
    
    
    p.def("get_id_instrument_curvAbs_table", []()
    {
        msg_error("SofaPython3::Input") << "Implementation is missing";
        return std::pair<int, int>{1,5};
    },sofapython3::doc::controller::InterventionalRadiologyController::InterventionalRadiologyControllerClass);

    py::module InterventionalRadiologyControllerBind = m.def_submodule("InterventionalRadiologyControllerBind");
    
    
    p.def("fixFirstNodesWithUntil",
        [](InterventionalRadiologyController<Rigid3Types>* self, unsigned int firstSimulatedNode ){
            self->fixFirstNodesWithUntil(firstSimulatedNode);
        });

    /// TODO: fill the docstring!
    //     InterventionalRadiologyControllerBind.doc() = R"doc(
    //            InterventionalRadiologyController
    //            -----------------------
    // 
    //            Advanced timer, meant to gather precise statistics for results in published papers.
    //            Not so advanced for now, but it will be...
    //        )doc";
    
    


}

} // namespace sofapython3
