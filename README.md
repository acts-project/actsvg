# actsvg

An SVG based C++17 plotting library for ACTS detectors, surfaces and relations.

This library has itself no dependency, simply for unittesting it relies on `googletest`. For the applications and examples, additionally `Boost` is used.

## Core module

This module has the entire definition and plotting functionality. 

## Meta module

This module acts as a translation layer between the caller libraries (etc. `ACTS`, `detray`) and core library.
It allows to create `proto` objects for detectors that can then be used for plotting using the `Core` library.

## Sample SVGs that can be produced

<table>
<tr>
<td width=200><img src="https://github.com/acts-project/actsvg/blob/main/docs/svg/odd_pixel_barrel_xy.svg" width=200></td>
<td width=200><img src="https://github.com/acts-project/actsvg/blob/main/docs/svg/odd_pixel_ec_xy.svg" width=200></td>
<td width=200><img src="https://github.com/acts-project/actsvg/blob/main/docs/svg/odd_pixel_ec_grid_xy.svg" width=200></td>
</tr>
<tr>
<td width=200><img src="https://github.com/acts-project/actsvg/blob/main/docs/svg/basic_rectangle.svg" width=200></td>
<td width=200><img src="https://github.com/acts-project/actsvg/blob/main/docs/svg/basic_trapezoid.svg" width=200></td>
<td width=200></td>
</tr>
</table>
