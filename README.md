# q2tools-220
Q2 compiler tools featuring ability to bsp v220 map format, automatic phong and soft spotlights, and other enhancements.
Includes qbsp3, qvis3, qrad3, and qdata.

Based on modifications supporting the v220 map format by XaeroX and DarkEssence distributed with the J.A.C.K. level editor.

The code is 'alpha' - not much testing or feedback at this point.

Fixes and enhancements:
*   AA tools (Alien Arena)- File path determination asumptions:
    *   moddir is parent of whatever directory contains the .map/.bsp
    *   gamedir is parent of moddir
    *   qdir is parent of gamedir

*   AA tools
    *   Use any existing TGA replacement textures for radiosity
    *   Spotlight center to surface point attenuation

*   DWH- Fix for scaled textures using an origin brush.

*   GDD tools (Geoffrey DeWan)
    *   Load file from PAK
    *   Split microbrushes
    *   Radiosity texture checking
    *   Update thread handling

*   jit (Jitspoe)
    *   Fix microbrush deformation.
    *   Re-enable PLANESIDE_EPSILON

*   KDT- qdata LWO support
*   kmbsp3 (Knightmare bsp tool)
    *   Caulk

*   MaxEd
    *   Texinfo overflow check
    *   VS compiler fixes

*   qbism
    * edge lighting fix 

*   quemap
    *   Face extents

*   VHLT (vluzacn)
    *   Automatic phong smoothing
    *   Add face for vertext normal
