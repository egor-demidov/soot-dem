//
// Created by egor on 7/19/24.
//

#include "preview_renderer.h"

Renderer::Renderer(std::vector<Eigen::Vector3d> const & x, double r_part) : r_part{r_part} {
    colors = vtkNew<vtkNamedColors>();

    monomer_representations.reserve(x.size());

    std::for_each(x.begin(), x.end(), [this](auto const & x [[]]) {
        vtkNew<vtkSphereSource> sphere;
        sphere->SetRadius(1.0);

        vtkNew<vtkOpenGLPolyDataMapper> sphere_mapper;
        sphere_mapper->SetInputConnection(sphere->GetOutputPort());

        vtkNew<vtkOpenGLActor> sphere_actor;
        sphere_actor->SetMapper(sphere_mapper);

        this->monomer_representations.emplace_back(sphere, sphere_mapper, sphere_actor);
    });

    renderer = vtkNew<vtkOpenGLRenderer>();
    std::for_each(monomer_representations.begin(), monomer_representations.end(), [this](auto monomer) {
        auto [source, mapper, actor] = monomer;
        this->renderer->AddActor(actor);
    });

    window = vtkNew<vtkRenderWindow>();
    window->AddRenderer(renderer);
    window->SetSize(800, 800);

    dump_number_label = vtkNew<vtkTextActor>();
    dump_number_label->SetPosition(10, 10);
    dump_number_label->GetTextProperty()->SetFontSize(20);
    renderer->AddActor2D(dump_number_label);

    update_preview(x, 0);
}

Renderer::~Renderer() = default;

void Renderer::update_preview(std::vector<Eigen::Vector3d> const & x, size_t dump_number) {
    for (size_t i = 0; i < x.size(); i ++) {
        std::get<2>(monomer_representations[i])->SetPosition(x[i][0] / r_part, x[i][1] / r_part, x[i][2] / r_part);
    }
    std::string label_content = "Dump #" + std::to_string(dump_number);
    dump_number_label->SetInput(label_content.c_str());

    window->Render();
}
