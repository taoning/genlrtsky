# genlrtsky
A command-line tool that uses libradtran to generate spectral sky description for Radiance ray-tracing engine.

The program outputs a sky defined by a Radiance hyperspectra (.hsr) fisheye image with a spectral range of 380-780nm at 20nm interval.

## Requirements
1. [libradtran](https://www.libradtran.org/doku.php?id=start)
2. [Optical properties of clouds and OPAC aerosol in netcdf format.](https://www.libradtran.org/doku.php?id=download)

## Usage
The command-line interface of `genlrtsky` follows the patterns of other
sky generators in Radiance, such as `gensky` and `gendaylit`.
The interface is four positional arguments followed by optional arguments:

Usage: genlrtsky month day hour minute [options]

Options:
  -a float
        Latitude, north positive (default 37)
  -c float
        Cloud cover 0-1
  -d float
        Aerosol optical depth
  -g float
        Albedo (default 0.2)
  -m int
        Standard meridian, west positive (default 120)
  -o float
        Longitude, west positive (default 122)
  -p string
        Output file prefix (default "default_")
  -quiet
        Quiet mode
  -s string
        Standard aerosol profile name (default "ca")
  -version
        Print version
  -y int
        Year (default 2003)

Positional arguments:
  month  Specify the month (1-12)
  day    Specify the day (1-31)
  hour   Specify the hour (0-23)
  minute Specify the minute (0-59)

### Aerosol
We uses the standard aerosol profiles provided by the OPAC packages.
To use one of the profile, include flag `-s`, followed by the shorthanded names, as seen in parentheses below:

* continental_average (ca)
* continental_clean (cc)
* continental_polluted (cp)
* maritime_clear (mc)
* maritime_tropical (mt)
* maritime_polluted (mp)
* antarctic (a)
* urban (u)
* desert (d)
* desert_spheroids (ds)

There is also a way to scale the aerosol profile to single broadband aerosol optical depth value.
For example, we can include `-d 0.15` flag to scale the profile to 0.15. 


## Examples

Generate a winter solstice noon clear sky with `continental_average` standard aerosol profile
`genlrtsky 12 21 12 0 -a 37 -o 122 -m 120 -s cc`

