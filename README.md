# actsvg

An SVG based C++17 plotting library for ACTS detectors, surfaces and relations.

This library has itself no dependency, simply for unittesting it relies on `googletest`.
## Core module

This module has the entire definition and plotting functionality. 

## Meta module

This module acts as a translation layer between the caller libraries (etc. `ACTS`, `detray`) and core library.
It allows to create `proto` objects for detectors that can then be used for plotting using the `core` module.
This module does also contain a dedicated `display` area that prepares some standard sheets, a view examples can be placed below.

## Sample SVGs that can be produced

A view stand-alone files that are produced by the `googletest` based unit testing.

<table>
<tr>
<td width=200><img src="https://github.com/acts-project/actsvg/blob/main/docs/svg/odd_pixel_barrel_xy.svg" width=200></td>
<td width=200><img src="https://github.com/acts-project/actsvg/blob/main/docs/svg/odd_pixel_endcap_xy.svg" width=200></td>
<td width=200><img src="https://github.com/acts-project/actsvg/blob/main/docs/svg/odd_pixel_endcap_grid_xy.svg" width=200></td>
</tr>
<tr>
<td width=200><img src="https://github.com/acts-project/actsvg/blob/main/docs/svg/basic_rectangle.svg" width=200></td>
<td width=200><img src="https://github.com/acts-project/actsvg/blob/main/docs/svg/basic_trapezoid.svg" width=200></td>
<td width=200></td>
</tr>
</table>

An example for an endcap sheet:

<img src="https://github.com/acts-project/actsvg/blob/main/docs/svg/odd_endcap_sheet_module_info.svg" width=600/>

An example for a barrel sheet (module info display):

<img src="https://github.com/acts-project/actsvg/blob/main/docs/svg/odd_barrel_sheet_module_info.svg" width=600/>