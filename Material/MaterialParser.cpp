////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2017-2018 Theodore Chang
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
////////////////////////////////////////////////////////////////////////////////

#include <Domain/DomainBase.h>
#include <Domain/ExternalModule.h>
#include <Material/Material>
#include <Toolbox/utility.h>

int create_new_material(const shared_ptr<DomainBase>& domain, istringstream& command) {
    string material_id;
    if(!get_input(command, material_id)) {
        suanpan_info("create_new_material() needs a tag.\n");
        return 0;
    }

    unique_ptr<Material> new_material = nullptr;

    if(is_equal(material_id, "Bilinear1D"))
        new_bilinear1d(new_material, command);
    else if(is_equal(material_id, "Bilinear2D"))
        new_bilinear2d(new_material, command);
    else if(is_equal(material_id, "Bilinear3D"))
        new_bilinear3d(new_material, command);
    else if(is_equal(material_id, "BilinearElastic1D"))
        new_bilinearelastic1d(new_material, command);
    else if(is_equal(material_id, "CDP"))
        new_cdp(new_material, command);
    else if(is_equal(material_id, "CDPM2"))
        new_cdpm2(new_material, command);
    else if(is_equal(material_id, "CDPPS"))
        new_cdpps(new_material, command);
    else if(is_equal(material_id, "Concrete01"))
        new_concrete01(new_material, command);
    else if(is_equal(material_id, "Concrete02"))
        new_concrete02(new_material, command);
    else if(is_equal(material_id, "Concrete21"))
        new_concrete21(new_material, command);
    else if(is_equal(material_id, "Concrete22"))
        new_concrete22(new_material, command);
    else if(is_equal(material_id, "Elastic1D"))
        new_elastic1d(new_material, command);
    else if(is_equal(material_id, "Elastic2D"))
        new_elastic2d(new_material, command);
    else if(is_equal(material_id, "Elastic3D"))
        new_elastic3d(new_material, command);
    else if(is_equal(material_id, "Gap01"))
        new_gap01(new_material, command);
    else if(is_equal(material_id, "Maxwell"))
        new_maxwell(new_material, command);
    else if(is_equal(material_id, "MooneyRivlin"))
        new_mooneyrivlin(new_material, command);
    else if(is_equal(material_id, "MPF"))
        new_mpf(new_material, command);
    else if(is_equal(material_id, "PlaneStrain"))
        new_planestrain(new_material, command);
    else if(is_equal(material_id, "PlaneStress"))
        new_planestress(new_material, command);
    else if(is_equal(material_id, "RambergOsgood"))
        new_rambergosgood(new_material, command);
    else if(is_equal(material_id, "RC01"))
        new_rc01(new_material, command);
    else if(is_equal(material_id, "Rebar2D"))
        new_rebar2d(new_material, command);
    else if(is_equal(material_id, "Viscosity1D"))
        new_viscosity1d(new_material, command);
    else {
        // check if the library is already loaded
        auto code = 0;
        for(const auto& I : domain->get_external_module_pool())
            if(I->library_name == material_id) {
                code = 1;
                break;
            }

        // not loaded then try load it
        if(code == 0 && domain->insert(make_shared<ExternalModule>(material_id))) code = 1;

        // if loaded find corresponding function
        if(code == 1)
            for(const auto& I : domain->get_external_module_pool()) {
                if(I->locate_module(material_id)) I->new_object(new_material, command);
                if(new_material != nullptr) break;
            }
    }

    if(new_material == nullptr || !domain->insert(move(new_material))) suanpan_debug("create_new_material() fails to insert new material.\n");

    return 0;
}

void new_bilinear1d(unique_ptr<Material>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_bilinear1d() requires a valid tag.\n");
        return;
    }

    double elastic_modulus;
    if(!get_input(command, elastic_modulus)) {
        suanpan_error("new_bilinear1d() requires a valid elastic modulus.\n");
        return;
    }

    double yield_stress;
    if(!get_input(command, yield_stress)) {
        suanpan_error("new_bilinear1d() requires a valid yield stress.\n");
        return;
    }

    auto hardening_ratio = 0.;
    if(!command.eof()) {
        if(!get_input(command, hardening_ratio)) {
            suanpan_error("new_bilinear1d() requires a valid hardening ratio.\n");
            return;
        }
    } else
        suanpan_debug("new_bilinear1d() assumes zero hardening ratio.\n");

    auto beta = 0.;
    if(!command.eof()) {
        if(!get_input(command, beta)) {
            suanpan_error("new_bilinear1d() requires a valid beta.\n");
            return;
        }
    } else
        suanpan_debug("new_bilinear1d() assumes isotropic hardening.\n");
    if(beta > 1.)
        beta = 1.;
    else if(beta < 0.)
        beta = 0.;

    auto density = 0.;
    if(!command.eof()) {
        if(!get_input(command, density)) {
            suanpan_error("new_bilinear1d() requires a valid density.\n");
            return;
        }
    } else
        suanpan_debug("new_bilinear1d() assumes zero density.\n");

    return_obj = make_unique<Bilinear1D>(tag, elastic_modulus, yield_stress, hardening_ratio, beta, density);
}

void new_bilinear2d(unique_ptr<Material>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_bilinear2d() requires a valid tag.\n");
        return;
    }

    double elastic_modulus;
    if(!get_input(command, elastic_modulus)) {
        suanpan_error("new_bilinear2d() requires a valid elastic modulus.\n");
        return;
    }

    double poissons_ratio;
    if(!get_input(command, poissons_ratio)) {
        suanpan_error("new_bilinear2d() requires a valid poissons ratio.\n");
        return;
    }

    double yield_stress;
    if(!get_input(command, yield_stress)) {
        suanpan_error("new_bilinear2d() requires a valid yield stress.\n");
        return;
    }

    auto hardening_ratio = 0.;
    if(!command.eof()) {
        if(!get_input(command, hardening_ratio)) {
            suanpan_error("new_bilinear2d() requires a valid hardening ratio.\n");
            return;
        }
    } else
        suanpan_debug("new_bilinear2d() assumes zero hardening ratio.\n");

    auto beta = 0.;
    if(!command.eof()) {
        if(!get_input(command, beta)) {
            suanpan_error("new_bilinear2d() requires a valid beta.\n");
            return;
        }
    } else
        suanpan_debug("new_bilinear2d() assumes isotropic hardening.\n");

    unsigned material_type = 0;
    if(!command.eof()) {
        if(!get_input(command, material_type)) {
            suanpan_error("new_bilinear2d() requires a valid material type.\n");
            return;
        }
    }

    auto density = 0.;
    if(!command.eof()) {
        if(!get_input(command, density)) {
            suanpan_error("new_bilinear2d() requires a valid density.\n");
            return;
        }
    } else
        suanpan_debug("new_bilinear2d() assumes zero density.\n");

    return_obj = make_unique<Bilinear2D>(tag, elastic_modulus, poissons_ratio, yield_stress, hardening_ratio, beta, material_type == 0 ? PlaneType::S : PlaneType::E, density);
}

void new_bilinear3d(unique_ptr<Material>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_bilinear3d() requires a valid tag.\n");
        return;
    }

    double elastic_modulus;
    if(!get_input(command, elastic_modulus)) {
        suanpan_error("new_bilinear3d() requires a valid elastic modulus.\n");
        return;
    }

    double poissons_ratio;
    if(!get_input(command, poissons_ratio)) {
        suanpan_error("new_bilinear3d() requires a valid poissons ratio.\n");
        return;
    }

    double yield_stress;
    if(!get_input(command, yield_stress)) {
        suanpan_error("new_bilinear3d() requires a valid yield stress.\n");
        return;
    }

    auto hardening_ratio = 0.;
    if(!command.eof()) {
        if(!get_input(command, hardening_ratio)) {
            suanpan_error("new_bilinear3d() requires a valid hardening ratio.\n");
            return;
        }
    } else
        suanpan_debug("new_bilinear3d() assumes zero hardening ratio.\n");

    auto beta = 0.;
    if(!command.eof()) {
        if(!get_input(command, beta)) {
            suanpan_error("new_bilinear3d() requires a valid beta.\n");
            return;
        }
    } else
        suanpan_debug("new_bilinear3d() assumes isotropic hardening.\n");

    auto density = 0.;
    if(!command.eof()) {
        if(!get_input(command, density)) {
            suanpan_error("new_bilinear3d() requires a valid density.\n");
            return;
        }
    } else
        suanpan_debug("new_bilinear3d() assumes zero density.\n");

    return_obj = make_unique<Bilinear3D>(tag, elastic_modulus, poissons_ratio, yield_stress, hardening_ratio, beta, density);
}

void new_bilinearelastic1d(unique_ptr<Material>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_bilinearelastic1d() requires a valid tag.\n");
        return;
    }

    double elastic_modulus;
    if(!get_input(command, elastic_modulus)) {
        suanpan_error("new_bilinearelastic1d() requires a valid elastic modulus.\n");
        return;
    }

    double yield_stress;
    if(!get_input(command, yield_stress)) {
        suanpan_error("new_bilinearelastic1d() requires a valid yield stress.\n");
        return;
    }

    auto hardening_ratio = 0.;
    if(!command.eof()) {
        if(!get_input(command, hardening_ratio)) {
            suanpan_error("new_bilinearelastic1d() requires a valid hardening ratio.\n");
            return;
        }
    } else
        suanpan_debug("new_bilinearelastic1d() assumes zero hardening ratio.\n");

    auto density = 0.;
    if(!command.eof()) {
        if(!get_input(command, density)) {
            suanpan_error("new_bilinearelastic1d() requires a valid density.\n");
            return;
        }
    } else
        suanpan_debug("new_bilinear1d() assumes zero density.\n");

    return_obj = make_unique<BilinearElastic1D>(tag, elastic_modulus, yield_stress, hardening_ratio, density);
}

void new_cdp(unique_ptr<Material>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_cdp() requires a valid tag.\n");
        return;
    }
    vec para_pool(std::initializer_list<double>{ 3E4, .2, 3., 30., 5E-4, 5E-2, .2, 2., .5, .65, .2, 1.16, .5, 2400E-12 });

    auto idx = 0;
    double para;
    while(!command.eof() && idx < 14)
        if(get_input(command, para)) para_pool(idx++) = para;

    return_obj = make_unique<CDP>(tag, para_pool(0), para_pool(1), para_pool(2), para_pool(3), para_pool(4), para_pool(5), para_pool(6), para_pool(7), para_pool(8), para_pool(9), para_pool(10), para_pool(11), para_pool(12), para_pool(13));
}

void new_cdpm2(unique_ptr<Material>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_cdp() requires a valid tag.\n");
        return;
    }
    return_obj = make_unique<CDPM2>(tag);
}

void new_cdpps(unique_ptr<Material>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_cdpps() requires a valid tag.\n");
        return;
    }
    vec para_pool(std::initializer_list<double>{ 3E4, .2, 3., 30., 5E-4, 5E-2, .2, 2., .5, .65, .2, 1.16, .5, 2400E-12 });

    auto idx = 0;
    double para;
    while(!command.eof() && idx < 14)
        if(get_input(command, para)) para_pool(idx++) = para;

    return_obj = make_unique<CDPPS>(tag, para_pool(0), para_pool(1), para_pool(2), para_pool(3), para_pool(4), para_pool(5), para_pool(6), para_pool(7), para_pool(8), para_pool(9), para_pool(10), para_pool(11), para_pool(12), para_pool(13));
}

void new_concrete01(unique_ptr<Material>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_concrete01() requires a valid tag.\n");
        return;
    }

    double peak_c_stress;
    if(!get_input(command, peak_c_stress)) {
        suanpan_error("new_concrete01() requires a valid compression stress.\n");
        return;
    }

    string backbone_type = "tsai";
    if(!command.eof() && !get_input(command, backbone_type)) {
        suanpan_error("new_concrete01() requires a valid backbone type.\n");
        return;
    }

    string center_oriented = "false";
    if(!command.eof() && !get_input(command, center_oriented)) {
        suanpan_error("new_concrete01() requires a valid center oriented switch.\n");
        return;
    }

    string no_tension = "false";
    if(!command.eof() && !get_input(command, no_tension)) {
        suanpan_error("new_concrete01() requires a valid tension behaviour switch.\n");
        return;
    }

    string tension_type = "linear";
    if(!command.eof() && !get_input(command, tension_type)) {
        suanpan_error("new_concrete01() requires a valid tension softening type.\n");
        return;
    }

    auto gf = 5E-3;
    if(!command.eof() && !get_input(command, gf)) {
        suanpan_error("new_concrete01() requires a valid fracture energy.\n");
        return;
    }

    auto density = 0.;
    if(!command.eof()) {
        if(!get_input(command, density)) {
            suanpan_error("new_concrete01() requires a valid density.\n");
            return;
        }
    } else
        suanpan_extra_debug("new_concrete01() assumes zero density.\n");

    BackboneType type;
    if(is_equal(backbone_type, "THORENFELDT"))
        type = BackboneType::THORENFELDT;
    else if(is_equal(backbone_type, "POPOVICS"))
        type = BackboneType::POPOVICS;
    else if(is_equal(backbone_type, "TSAI"))
        type = BackboneType::TSAI;
    else if(is_equal(backbone_type, "KPS"))
        type = BackboneType::KPS;
    else {
        suanpan_error("new_concrete01() cannot identify backbone type.\n");
        return;
    }

    TensionType typea;
    if(is_equal(tension_type, "EXP"))
        typea = TensionType::EXPONENTIAL;
    else if(is_equal(tension_type, "EXPONENTIAL"))
        typea = TensionType::EXPONENTIAL;
    else if(is_equal(tension_type, "LINEAR"))
        typea = TensionType::LINEAR;
    else {
        suanpan_error("new_concrete01() cannot identify tension softening type.\n");
        return;
    }

    return_obj = make_unique<Concrete01>(tag, peak_c_stress, type, is_true(center_oriented), is_true(no_tension), typea, gf, density);
}

void new_concrete02(unique_ptr<Material>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_concrete01() requires a valid tag.\n");
        return;
    }

    double peak_c_strain;
    if(!get_input(command, peak_c_strain)) {
        suanpan_error("new_concrete01() requires a valid compression strain.\n");
        return;
    }

    double peak_c_stress;
    if(!get_input(command, peak_c_stress)) {
        suanpan_error("new_concrete01() requires a valid compression stress.\n");
        return;
    }

    auto density = 0.;
    if(!command.eof()) {
        if(!get_input(command, density)) {
            suanpan_error("new_concrete01() requires a valid density.\n");
            return;
        }
    } else
        suanpan_extra_debug("new_concrete01() assumes zero density.\n");

    return_obj = make_unique<Concrete02>(tag, peak_c_stress, density);
}

void new_concrete21(unique_ptr<Material>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_concrete21() requires a valid tag.\n");
        return;
    }

    double peak_c_stress;
    if(!get_input(command, peak_c_stress)) {
        suanpan_error("new_concrete21() requires a valid compression stress.\n");
        return;
    }

    string backbone_type;
    if(!get_input(command, backbone_type)) {
        suanpan_error("new_concrete21() requires a valid backbone type.\n");
        return;
    }

    string center_oriented = "false";
    if(!command.eof() && !get_input(command, center_oriented)) {
        suanpan_error("new_concrete21() requires a valid center oriented switch.\n");
        return;
    }

    string tension_type = "linear";
    if(!command.eof() && !get_input(command, tension_type)) {
        suanpan_error("new_concrete21() requires a valid tension softening type.\n");
        return;
    }

    auto gf = 1E-2;
    if(!command.eof() && !get_input(command, gf)) {
        suanpan_error("new_concrete21() requires a valid fracture energy.\n");
        return;
    }

    string poisson = "true";
    if(!command.eof() && !get_input(command, poisson)) {
        suanpan_error("new_concrete21() requires a valid poisson switch.\n");
        return;
    }

    string degradation = "false";
    if(!command.eof() && !get_input(command, degradation)) {
        suanpan_error("new_concrete21() requires a valid degradation switch.\n");
        return;
    }

    auto density = 0.;
    if(!command.eof()) {
        if(!get_input(command, density)) {
            suanpan_error("new_concrete21() requires a valid density.\n");
            return;
        }
    } else
        suanpan_extra_debug("new_concrete21() assumes zero density.\n");

    BackboneType type;
    if(is_equal(backbone_type, "THORENFELDT"))
        type = BackboneType::THORENFELDT;
    else if(is_equal(backbone_type, "POPOVICS"))
        type = BackboneType::POPOVICS;
    else if(is_equal(backbone_type, "TSAI"))
        type = BackboneType::TSAI;
    else if(is_equal(backbone_type, "KPS"))
        type = BackboneType::KPS;
    else {
        suanpan_error("new_concrete21() cannot identify backbone type.\n");
        return;
    }

    TensionType typea;
    if(is_equal(tension_type, "EXP"))
        typea = TensionType::EXPONENTIAL;
    else if(is_equal(tension_type, "EXPONENTIAL"))
        typea = TensionType::EXPONENTIAL;
    else if(is_equal(tension_type, "LINEAR"))
        typea = TensionType::LINEAR;
    else {
        suanpan_error("new_concrete01() cannot identify tension softening type.\n");
        return;
    }

    return_obj = make_unique<Concrete21>(tag, peak_c_stress, type, is_true(center_oriented), typea, gf, is_true(poisson), is_true(degradation), density);
}

void new_concrete22(unique_ptr<Material>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_concrete22() requires a valid tag.\n");
        return;
    }

    double peak_c_stress;
    if(!get_input(command, peak_c_stress)) {
        suanpan_error("new_concrete22() requires a valid compression stress.\n");
        return;
    }

    string backbone_type;
    if(!get_input(command, backbone_type)) {
        suanpan_error("new_concrete22() requires a valid backbone type.\n");
        return;
    }

    auto shear_retention = .2;
    if(!command.eof() && (!get_input(command, shear_retention) || shear_retention < 0. || shear_retention > 1.)) {
        suanpan_error("new_concrete22() requires a valid shear retention factor.\n");
        return;
    }

    string center_oriented = "false";
    if(!command.eof() && !get_input(command, center_oriented)) {
        suanpan_error("new_concrete22() requires a valid center oriented switch.\n");
        return;
    }

    string tension_type = "linear";
    if(!command.eof() && !get_input(command, tension_type)) {
        suanpan_error("new_concrete21() requires a valid tension softening type.\n");
        return;
    }

    auto gf = 1E-2;
    if(!command.eof() && !get_input(command, gf)) {
        suanpan_error("new_concrete22() requires a valid fracture energy.\n");
        return;
    }

    string poisson = "true";
    if(!command.eof() && !get_input(command, poisson)) {
        suanpan_error("new_concrete22() requires a valid poisson switch.\n");
        return;
    }

    string degradation = "false";
    if(!command.eof() && !get_input(command, degradation)) {
        suanpan_error("new_concrete22() requires a valid degradation switch.\n");
        return;
    }

    auto density = 0.;
    if(!command.eof()) {
        if(!get_input(command, density)) {
            suanpan_error("new_concrete22() requires a valid density.\n");
            return;
        }
    } else
        suanpan_extra_debug("new_concrete22() assumes zero density.\n");

    BackboneType type;
    if(is_equal(backbone_type, "THORENFELDT"))
        type = BackboneType::THORENFELDT;
    else if(is_equal(backbone_type, "POPOVICS"))
        type = BackboneType::POPOVICS;
    else if(is_equal(backbone_type, "TSAI"))
        type = BackboneType::TSAI;
    else if(is_equal(backbone_type, "KPS"))
        type = BackboneType::KPS;
    else {
        suanpan_error("new_concrete22() cannot identify backbone type.\n");
        return;
    }

    TensionType typea;
    if(is_equal(tension_type, "EXP"))
        typea = TensionType::EXPONENTIAL;
    else if(is_equal(tension_type, "EXPONENTIAL"))
        typea = TensionType::EXPONENTIAL;
    else if(is_equal(tension_type, "LINEAR"))
        typea = TensionType::LINEAR;
    else {
        suanpan_error("new_concrete01() cannot identify tension softening type.\n");
        return;
    }

    return_obj = make_unique<Concrete22>(tag, peak_c_stress, type, shear_retention, is_true(center_oriented), typea, gf, is_true(poisson), is_true(degradation), density);
}

void new_elastic1d(unique_ptr<Material>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_elastic1d() requires a valid tag.\n");
        return;
    }

    double elastic_modulus;
    if(!get_input(command, elastic_modulus)) {
        suanpan_error("new_elastic1d() requires a valid elastic modulus.\n");
        return;
    }

    auto density = 0.;
    if(!command.eof()) {
        if(!get_input(command, density)) {
            suanpan_error("new_elastic1d() requires a valid density.\n");
            return;
        }
    } else
        suanpan_debug("new_elastic1d() assumes zero density.\n");

    return_obj = make_unique<Elastic1D>(tag, elastic_modulus, density);
}

void new_elastic2d(unique_ptr<Material>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_elastic2d() requires a valid tag.\n");
        return;
    }

    double elastic_modulus;
    if(!get_input(command, elastic_modulus)) {
        suanpan_error("new_elastic2d() requires a valid elastic modulus.\n");
        return;
    }

    double poissons_ratio;
    if(!get_input(command, poissons_ratio)) {
        suanpan_error("new_elastic2d() requires a valid poissons ratio.\n");
        return;
    }

    auto density = 0.;
    if(!command.eof()) {
        if(!get_input(command, density)) {
            suanpan_error("new_elastic2d() requires a valid density.\n");
            return;
        }
    } else
        suanpan_debug("new_elastic2d() assumes zero density.\n");

    auto material_type = 0;
    if(!command.eof()) {
        if(!get_input(command, material_type)) {
            suanpan_error("new_elastic2d() requires a valid material type.\n");
            return;
        }
    } else
        suanpan_debug("new_elastic2d() assumes plane stress.\n");

    return_obj = make_unique<Elastic2D>(tag, elastic_modulus, poissons_ratio, density, material_type == 0 ? PlaneType::S : PlaneType::E);
}

void new_elastic3d(unique_ptr<Material>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_elastic3d() requires a valid tag.\n");
        return;
    }

    double elastic_modulus;
    if(!get_input(command, elastic_modulus)) {
        suanpan_error("new_elastic3d() requires a valid elastic modulus.\n");
        return;
    }

    double poissons_ratio;
    if(!get_input(command, poissons_ratio)) {
        suanpan_error("new_elastic3d() requires a valid poissons ratio.\n");
        return;
    }

    auto density = 0.;
    if(!command.eof()) {
        if(!get_input(command, density)) {
            suanpan_error("new_elastic3d() requires a valid density.\n");
            return;
        }
    } else
        suanpan_debug("new_elastic3d() assumes zero density.\n");

    return_obj = make_unique<Elastic3D>(tag, elastic_modulus, poissons_ratio, density);
}

void new_gap01(unique_ptr<Material>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_gap01() requires a valid tag.\n");
        return;
    }

    double elastic_modulus;
    if(!get_input(command, elastic_modulus)) {
        suanpan_error("new_gap01() requires a valid elastic modulus.\n");
        return;
    }

    double yield_stress;
    if(!get_input(command, yield_stress)) {
        suanpan_error("new_gap01() requires a valid yield stress.\n");
        return;
    }

    auto gap_strain = 0.;
    if(!command.eof() && !get_input(command, gap_strain)) {
        suanpan_error("new_gap01() requires a valid hardening ratio.\n");
        return;
    }

    auto density = 0.;
    if(!command.eof() && !get_input(command, density)) {
        suanpan_error("new_gap01() requires a valid density.\n");
        return;
    }

    return_obj = make_unique<Gap01>(tag, elastic_modulus, yield_stress, gap_strain, density);
}

void new_maxwell(unique_ptr<Material>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_maxwell() requires a valid tag.\n");
        return;
    }

    double elastic_modulus;
    if(!get_input(command, elastic_modulus)) {
        suanpan_error("new_maxwell() requires a valid elastic modulus.\n");
        return;
    }

    double alpha;
    if(!get_input(command, alpha)) {
        suanpan_error("new_maxwell() requires a valid alpha.\n");
        return;
    }

    double damping_postive;
    if(!get_input(command, damping_postive)) {
        suanpan_error("new_maxwell() requires a valid damping coefficient.\n");
        return;
    }

    return_obj = make_unique<Maxwell>(tag, elastic_modulus, alpha, damping_postive);
}

void new_mooneyrivlin(unique_ptr<Material>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_mooneyrivlin() requires a valid tag.\n");
        return;
    }

    double bulk_modulus;
    if(!get_input(command, bulk_modulus)) {
        suanpan_error("new_mooneyrivlin() requires a valid bulk modulus.\n");
        return;
    }

    double a10;
    if(!get_input(command, a10)) {
        suanpan_error("new_mooneyrivlin() requires a valid a10.\n");
        return;
    }

    double a01;
    if(!get_input(command, a01)) {
        suanpan_error("new_mooneyrivlin() requires a valid a01.\n");
        return;
    }

    auto density = 0.;
    if(!command.eof()) {
        if(!get_input(command, density)) {
            suanpan_error("new_mooneyrivlin() requires a valid density.\n");
            return;
        }
    } else
        suanpan_debug("new_mooneyrivlin() assumes zero density.\n");

    return_obj = make_unique<MooneyRivlin>(tag, bulk_modulus, a10, a01, density);
}

void new_mpf(unique_ptr<Material>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_mpf() requires a valid tag.\n");
        return;
    }

    double elastic_modulus;
    if(!get_input(command, elastic_modulus)) {
        suanpan_error("new_mpf() requires a valid elastic modulus.\n");
        return;
    }

    double yield_stress;
    if(!get_input(command, yield_stress)) {
        suanpan_error("new_mpf() requires a valid yield stress.\n");
        return;
    }

    auto hardening_ratio = .05;
    if(!command.eof() && !get_input(command, hardening_ratio)) {
        suanpan_error("new_mpf() requires a valid hardening ratio.\n");
        return;
    }

    auto R0 = 20.;
    if(!command.eof() && !get_input(command, R0)) {
        suanpan_error("new_mpf() requires a valid R0.\n");
        return;
    }

    auto A1 = 18.5;
    if(!command.eof() && !get_input(command, A1)) {
        suanpan_error("new_mpf() requires a valid A1.\n");
        return;
    }

    auto A2 = .15;
    if(!command.eof() && !get_input(command, A2)) {
        suanpan_error("new_mpf() requires a valid A2.\n");
        return;
    }

    auto A3 = .01;
    if(!command.eof() && !get_input(command, A3)) {
        suanpan_error("new_mpf() requires a valid A3.\n");
        return;
    }

    auto A4 = 7.;
    if(!command.eof() && !get_input(command, A4)) {
        suanpan_error("new_mpf() requires a valid A4.\n");
        return;
    }

    string iso = "false";
    if(!command.eof() && !get_input(command, iso)) {
        suanpan_error("new_mpf() requires a valid isotropic hardening switch.\n");
        return;
    }

    string con = "false";
    if(!command.eof() && !get_input(command, con)) {
        suanpan_error("new_mpf() requires a valid constant radius switch.\n");
        return;
    }

    auto density = 0.;
    if(!command.eof() && !get_input(command, density)) {
        suanpan_error("new_mpf() requires a valid density.\n");
        return;
    }

    return_obj = make_unique<MPF>(tag, elastic_modulus, yield_stress, hardening_ratio, R0, A1, A2, A3, A4, is_true(iso), is_true(con), density);
}

void new_planestrain(unique_ptr<Material>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_planestrain() requires a valid tag.\n");
        return;
    }

    unsigned full_tag;
    if(!get_input(command, full_tag)) {
        suanpan_error("new_planestrain() requires a valid reference material tag.\n");
        return;
    }

    return_obj = make_unique<PlaneStrain>(tag, full_tag);
}

void new_planestress(unique_ptr<Material>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_planestress() requires a valid tag.\n");
        return;
    }

    unsigned full_tag;
    if(!get_input(command, full_tag)) {
        suanpan_error("new_planestress() requires a valid reference material tag.\n");
        return;
    }

    auto max_iteration = 20;
    if(!command.eof() && !get_input(command, max_iteration)) {
        suanpan_error("new_planestress() requires a number for maximum iteration.\n");
        return;
    }

    auto tolerance = 1E-10;
    if(!command.eof() && !get_input(command, tolerance)) {
        suanpan_error("new_planestress() requires a valid tolerance.\n");
        return;
    }

    string use_matrix;
    if(!command.eof() && !get_input(command, use_matrix)) {
        suanpan_error("new_planestress() requires a valid flag to indicate if to use the matrix in iteration.\n");
        return;
    }

    return_obj = make_unique<PlaneStress>(tag, full_tag, max_iteration, tolerance, is_true(use_matrix));
}

void new_rambergosgood(unique_ptr<Material>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_rambergosgood() requires a valid tag.\n");
        return;
    }

    double elastic_modulus;
    if(!get_input(command, elastic_modulus)) {
        suanpan_error("new_rambergosgood() requires a valid elastic modulus.\n");
        return;
    }

    double yield_stress;
    if(!get_input(command, yield_stress)) {
        suanpan_error("new_rambergosgood() requires a valid yield stress.\n");
        return;
    }

    auto offset = 1;
    if(!command.eof() && !get_input(command, offset)) {
        suanpan_error("new_rambergosgood() requires a valid offset.\n");
        return;
    }

    auto n = 10.;
    if(!command.eof() && !get_input(command, n)) {
        suanpan_error("new_rambergosgood() requires a valid n.\n");
        return;
    }

    auto density = 0.;
    if(!command.eof() && !get_input(command, density)) {
        suanpan_error("new_rambergosgood() requires a valid density.\n");
        return;
    }

    return_obj = make_unique<RambergOsgood>(tag, elastic_modulus, yield_stress, offset, n, density);
}

void new_rc01(unique_ptr<Material>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_rc01() requires a valid tag.\n");
        return;
    }

    unsigned rebar_tag, concrete_tag;
    if(!get_input(command, rebar_tag)) {
        suanpan_error("new_rc01() requires a valid rebar tag.\n");
        return;
    }
    if(!get_input(command, concrete_tag)) {
        suanpan_error("new_rc01() requires a valid concrete tag.\n");
        return;
    }

    return_obj = make_unique<RC01>(tag, rebar_tag, concrete_tag);
}

void new_rebar2d(unique_ptr<Material>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_rebarlayer() requires a valid tag.\n");
        return;
    }

    unsigned major_tag, minor_tag;
    if(!get_input(command, major_tag)) {
        suanpan_error("new_rebarlayer() requires a valid material tag.\n");
        return;
    }
    if(!get_input(command, minor_tag)) {
        suanpan_error("new_rebarlayer() requires a valid material tag.\n");
        return;
    }

    double major_ratio, minor_ratio;
    if(!get_input(command, major_ratio)) {
        suanpan_error("new_rebarlayer() requires a valid reinforcement ratio.\n");
        return;
    }
    if(!get_input(command, minor_ratio)) {
        suanpan_error("new_rebarlayer() requires a valid reinforcement ratio.\n");
        return;
    }

    return_obj = make_unique<Rebar2D>(tag, major_tag, minor_tag, major_ratio, minor_ratio);
}

void new_viscosity1d(unique_ptr<Material>& return_obj, istringstream& command) {
    unsigned tag;
    if(!get_input(command, tag)) {
        suanpan_error("new_viscosity1d() requires a valid tag.\n");
        return;
    }

    double alpha;
    if(!get_input(command, alpha)) {
        suanpan_error("new_viscosity1d() requires a valid alpha.\n");
        return;
    }

    double damping_a;
    if(!get_input(command, damping_a)) {
        suanpan_error("new_viscosity1d() requires a valid damping coefficient for the first quadrant.\n");
        return;
    }

    auto damping_b = damping_a;
    if(!command.eof() && !get_input(command, damping_b)) {
        suanpan_error("new_viscosity1d() requires a valid damping coefficient for the second quadrant.\n");
        return;
    }

    auto damping_c = damping_a;
    if(!command.eof() && !get_input(command, damping_c)) {
        suanpan_error("new_viscosity1d() requires a valid damping coefficient for the third quadrant.\n");
        return;
    }

    auto damping_d = damping_a;
    if(!command.eof() && !get_input(command, damping_d)) {
        suanpan_error("new_viscosity1d() requires a valid damping coefficient for the fourth quadrant.\n");
        return;
    }

    auto gap_a = 1E4;
    if(!command.eof() && !get_input(command, gap_a)) {
        suanpan_error("new_viscosity1d() requires a valid gap size for strain axis.\n");
        return;
    }

    auto gap_b = 1E4;
    if(!command.eof() && !get_input(command, gap_b)) {
        suanpan_error("new_viscosity1d() requires a valid gap size for strain rate axis.\n");
        return;
    }

    return_obj = make_unique<Viscosity1D>(tag, alpha, damping_a, damping_b, damping_c, damping_d, gap_a, gap_b);
}

int test_material(const shared_ptr<DomainBase>& domain, istringstream& command) {
    unsigned material_tag;
    if(!get_input(command, material_tag)) {
        suanpan_error("test_material() needs a valid material tag.\n");
        return 0;
    }

    double step_size;
    if(!get_input(command, step_size)) {
        suanpan_error("test_material() needs a valid step size.\n");
        return 0;
    }

    vector<unsigned> load_step;
    int step;
    while(get_input(command, step)) load_step.push_back(step < 0 ? -step : step);

    auto& material_proto = domain->get_material(material_tag);

    if(material_proto == nullptr) return 0;

    auto result = material_tester(material_proto->get_copy(), load_step, step_size);

    result.save("RESULT.h5", hdf5_binary_trans);

    return 0;
}
