// This file is part of the actsvg packge.
//
// Copyright (C) 2023 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <pybind11/eval.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <memory>
#include <tuple>

#include "actsvg/core.hpp"
#include "actsvg/meta.hpp"
#include "utilities.hpp"

namespace py = pybind11;
using namespace pybind11::literals;

namespace actsvg {
namespace python {

using surface = proto::surface<point3_collection>;

namespace {

/// @brief  Helper method to view surfaces given a view type
/// @tparam view_type
///
/// @param ss the surfaces, @note they will be copied to apply the style
/// @param apply_style is the boolean weather the style should be applied or not
/// @param f is the fill object
/// @param str is the stroke object
///
/// @return the surfaces
template <typename view_type>
auto view_surfaces(const std::vector<surface>& ss, bool apply_style = true,
                   const style::fill& f = __s_fill,
                   const style::stroke& str = __s_stroke) {

    view_type vw;
    // Loop over the surfaces and view them
    std::vector<svg::object> surfaces;
    surfaces.reserve(ss.size());
    for (const auto& s : ss) {
        // Make a copy for setting new parameters
        if (apply_style) {
            auto sc = s;
            sc._fill = f;
            sc._stroke = str;
            surfaces.push_back(display::surface(sc._name, sc, vw));
            continue;
        }
        // Add the non-copied one
        surfaces.push_back(display::surface(s._name, s, vw));
    }
    return surfaces;
}

/// draw the surface as polygon
///
/// @param s the surface
/// @param f the fill
/// @param str the stroke
/// @param view the view string
auto surface_as_oriented_polygon(const surface& s, const style::fill& f,
                                 const style::stroke& str,
                                 const std::string& view) {

    auto sc = s;
    sc._fill = f;
    sc._stroke = str;
    if (view == "xy") {
        views::x_y xy;
        return display::oriented_polygon(sc._name, s, xy);
    }
    if (view == "zphi") {
        views::z_phi zphi;
        return display::oriented_polygon(sc._name, s, zphi);
    }
    throw std::runtime_error(
        "surface_as_oriented_polygon: No matching view type found");
    // Return a sterile object
    return svg::object{};
}

/// @brief cast helper
///
/// @param v the vertex that should be casted
/// @param val the cast value
///
/// @return a single range
scalar cast(const point3& v, const std::string& val) {
    scalar vc = std::numeric_limits<scalar>::max();
    // Specify the cast selection
    if (val == "x") {
        vc = v[0u];
    } else if (val == "y") {
        vc = v[1u];
    } else if (val == "z") {
        vc = v[2u];
    } else if (val == "r") {
        vc = std::sqrt(v[0u] * v[0u] + v[1u] * v[1u]);
    } else if (val == "phi") {
        vc = std::atan2(v[1u], v[2u]);
    }
    return vc;
}

}  // namespace

/// @brief  Adding the proto module to the
/// python bindings
///
/// @param ctx the python context
void add_display_module(context& ctx) {

    auto& m = ctx.get("main");
    auto d = m.def_submodule("display");

    {
        /// Create a default style tuple
        d.def("create_default_style", []() {
            return std::tuple<style::fill, style::stroke>{__s_fill, __s_stroke};
        });
    }

    {
        /// Create a default style tuple
        /// @param f_rgb the fill color
        /// @param f_opac the fill opacity
        /// @param str_rgb the stroke color
        /// @param str_opac the stroke opacity
        /// @param str_w the stroke witdth
        d.def("create_custom_style",
              [](const std::array<int, 3u>& f_rgb, scalar f_opac,
                 const std::array<int, 3u>& str_rgb, scalar str_opac,
                 scalar str_w) {
                  // The fill definition
                  style::fill f;
                  f._fc._rgb = f_rgb;
                  f._fc._opacity = f_opac;
                  // The stroke definition
                  style::stroke s;
                  s._sc._rgb = str_rgb;
                  s._sc._opacity = str_opac;
                  s._width = str_w;
                  // Return the tuple
                  return std::tuple<style::fill, style::stroke>{f, s};
              });
    }

    {  /// Select a vertex view range
        ///
        /// @param vs are the vertices to be shown
        /// @param restr are the restrictions
        ///
        /// @return the survivors
        d.def("select_view_range",
              [](const point3_collection& vs,
                 const std::vector<
                     std::tuple<std::string, std::array<scalar, 2u>>>& restr) {
                  point3_collection vs_restr;
                  for (const auto& v : vs) {
                      bool veto = false;
                      // Loop over the range restrictions and veto
                      for (const auto& r : restr) {
                          auto [val, range] = r;
                          scalar vc = cast(v, val);
                          if (vc <= range[0u] or vc >= range[1u]) {
                              veto = true;
                              break;
                          }
                      }
                      // If not veoted, take it
                      if (not veto) {
                          vs_restr.push_back(v);
                      }
                  }
                  return vs_restr;
              });
    }

    {
        /// Select surfaces according to their vertices in a specific z range
        ///
        /// @param ss is the input surface vector
        /// @param val is the value in which should be selected: x,y,z,r,phi,eta
        /// supported
        /// @param range is the range for the selection
        ///
        /// @return selected surfaces
        d.def("select_view_range",
              [](const std::vector<surface>& ss, const std::string& val,
                 const std::array<scalar, 2u>& range) {
                  std::vector<surface> sss;
                  for (const auto& s : ss) {
                      bool select = false;
                      for (const auto& v : s._vertices) {
                          // Cast to the view rangne
                          scalar vc = cast(v, val);
                          if (vc >= range[0u] and vc <= range[1u]) {
                              select = true;
                              break;
                          }
                      }
                      // Selection was triggered
                      if (select) {
                          sss.push_back(s);
                      }
                  }
                  return sss;
              });
    }

    {

        /// View a single surface as oriented polygon - surface style
        ///
        /// @param s is the surface
        /// @param f is the fill style
        /// @param str is the stroke style
        /// @param view is the view string (xy, zphi)
        ///
        /// @note the surfaces will be copied to apply the style selection
        ///
        /// @return an object (sterile if now matching view type is found)
        d.def("surface_as_oriented_polygon",
              [](const surface& s, const style::fill& f,
                 const style::stroke& str, const std::string& view) {
                  return surface_as_oriented_polygon(s, f, str, view);
              });

        /// View a single surface
        ///
        /// @param s is the surface
        /// @param f is the fill style
        /// @param str is the stroke style
        /// @param view is the view string (xy, zphi)
        ///
        /// @note the surfaces will be copied to apply the style selection
        ///
        /// @return an object (sterile if now matching view type is found)
        d.def("surface", [](const surface& s, const style::fill& f,
                            const style::stroke& str, const std::string& view) {
            auto sc = s;
            sc._fill = f;
            sc._stroke = str;
            if (view == "xy") {
                views::x_y xy;
                return display::surface(sc._name, s, xy);
            }
            if (view == "zphi") {
                views::z_phi zphi;
                return display::surface(sc._name, s, zphi);
            }
            // Return a sterile object
            return svg::object{};
        });

        /// View a single surface as oriented polygon - style applied
        ///
        /// @param s is the surface
        /// @param f is the fill style
        /// @param str is the stroke style
        /// @param view is the view string (xy, zphi)
        ///
        /// @note the surfaces will be copied to apply the style selection
        ///
        /// @return an object (sterile if now matching view type is found)
        d.def("surface_as_oriented_polygon",
              [](const surface& s, const style::fill& f,
                 const style::stroke& str, const std::string& view) {
                  return surface_as_oriented_polygon(s, f, str, view);
              });

        /// View surfaces -  with surface internal style
        ///
        /// @param ss the surfaces
        /// @param view the view string (xy, zphi)
        ///
        /// @return a list of objects
        d.def("surfaces",
              [](const std::vector<surface>& ss, const std::string& view) {
                  // xy view is chosen
                  if (view == "xy") {
                      return view_surfaces<views::x_y>(ss, false);
                  }
                  // zphi view is chosen
                  if (view == "zphi") {
                      return view_surfaces<views::z_phi>(ss, false);
                  }
                  return std::vector<svg::object>{};
              });

        /// View surfaces as oriented polygons -  with surface internal style
        ///
        /// @param ss is the surface collection
        /// @param view is the view string (xy, zphi)
        ///
        /// @note the surfaces will be copied to apply the style selection
        ///
        /// @return an object (sterile if now matching view type is found)
        d.def("surfaces_as_oriented_polygons",
              [](const std::vector<surface>& ss, const std::string& view) {
                  std::vector<svg::object> polygons;
                  for (const auto& s : ss) {
                      polygons.push_back(surface_as_oriented_polygon(
                          s, s._fill, s._stroke, view));
                  }
                  // Return a sterile object
                  return polygons;
              });

        /// View surfaces -  with style appied
        ///
        /// @param ss the surfaces
        /// @param f the fill style
        /// @param str the stroke style
        /// @param view the view string (xy, zphi)
        ///
        /// @return a list of objects
        d.def("surfaces",
              [](const std::vector<surface>& ss, const style::fill& f,
                 const style::stroke str, const std::string& view) {
                  // xy view is chosen
                  if (view == "xy") {
                      return view_surfaces<views::x_y>(ss, true, f, str);
                  }
                  // zphi view is chosen
                  if (view == "zphi") {
                      return view_surfaces<views::z_phi>(ss, true, f, str);
                  }
                  return std::vector<svg::object>{};
              });
    }

    /// View surfaces as oriented polygons - style applied
    ///
    /// @param ss the surfaces
    /// @param f is the fill style
    /// @param str is the stroke style
    /// @param view is the view string (xy, zphi)
    ///
    /// @note the surfaces will be copied to apply the style selection
    ///
    /// @return an object (sterile if now matching view type is found)
    d.def("surfaces_as_oriented_polygons", [](const std::vector<surface>& ss,
                                              const style::fill& f,
                                              const style::stroke& str,
                                              const std::string& view) {
        std::vector<svg::object> polygons;
        for (const auto& s : ss) {
            polygons.push_back(surface_as_oriented_polygon(s, f, str, view));
        }
        // Return a sterile object
        return polygons;
    });

    {
        /// View step_tracks -  with style appied
        ///
        /// @param sts the surfaces
        /// @param rsize is the radius size
        /// @param f the fill style
        /// @param str the stroke style
        /// @param view the view string (xy, zphi)
        ///
        /// @return a list of objects
        d.def("step_tracks", [](const std::vector<std::vector<point3>>& sts,
                                scalar rsize, const style::fill& f,
                                const style::stroke& str, std::string& view) {
            std::vector<std::vector<svg::object>> track_views;
            track_views.reserve(sts.size());
            // Loop and construct
            for (auto [it, st] : utils::enumerate(sts)) {
                std::vector<svg::object> track_view;
                for (auto [is, s] : utils::enumerate(st)) {
                    scalar v_x = 0.;
                    scalar v_y = 0.;
                    if (view == "xy") {
                        v_x = s[0u];
                        v_y = s[1u];
                    }
                    // Add to the track view
                    track_view.push_back(
                        draw::circle("track_" + std::to_string(it) + "_step" +
                                         std::to_string(is),
                                     {v_x, v_y}, rsize, f, str));
                }
                track_views.push_back(track_view);
                track_view.clear();
            }
            return track_views;
        });
    }

    {
        /// View surface material
        ///
        /// @param sm the surface material
        /// @param view the view string (xy, zphi)
        ///
        /// @return a list of objects
        d.def("surface_material", &display::surface_material);
    }

}  // namespace actsvg

}  // namespace python
}  // namespace actsvg
