//
// Created by egor on 7/19/24.
//

#ifndef SOOT_AFM_PREVIEW_RENDERER_H
#define SOOT_AFM_PREVIEW_RENDERER_H

#include <vector>

#include <Eigen/Eigen>

#include <vtkOpenGLActor.h>
#include <vtkBoxWidget.h>
#include <vtkCamera.h>
#include <vtkCommand.h>
#include <vtkSphereSource.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkOpenGLPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkOpenGLRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkOpenGLRenderer.h>
#include <vtkTransform.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>

class Renderer {
public:
    explicit Renderer(std::vector<Eigen::Vector3d> const & x, double r_part);
    ~Renderer();

    void update_preview(std::vector<Eigen::Vector3d> const & x, size_t dump_number);

private:
    typedef std::tuple<
        vtkSmartPointer<vtkOpenGLPolyDataMapper>,
        vtkSmartPointer<vtkOpenGLActor>> monomer_representation_t;

    const double r_part;
    vtkSmartPointer<vtkNamedColors> colors;
    vtkSmartPointer<vtkSphereSource> sphere_source;
    vtkSmartPointer<vtkOpenGLRenderer> renderer;
    vtkSmartPointer<vtkRenderWindow> window;
    vtkSmartPointer<vtkTextActor> dump_number_label;
    std::vector<monomer_representation_t> monomer_representations;
};

#endif //SOOT_AFM_PREVIEW_RENDERER_H
