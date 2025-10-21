```````
//==== SHARK: Soft contributions to Heavy quark production in ARbitrary Kinematics ====//
//                                     ,-                                              //
//                                   ,'::|                                             //
//                                  /::::|                                             //
//                                ,':::::\                                      _..    //
//             ____........-------...:::::\                                  ,-' /     //
//     _.--"""". . . .      .   .  .  .    `-._                           ,-' .;'      //
//    <. - :::::o......  ...   . . .. . .  .  .""--._                  ,-'. .;'	       //
//     `-._  ` `":`:`:`::||||:::::::::::::::::.:. .  ""--._ ,'|     ,-'.  .;'          //
//         """_=--       //    .. ````:`:`::::::::::.:.:.:. .`-`._-'.   .;'            //
//             ""--.__      (       \               ` ``:`:``:::: .   .;'              //
//                    "\""--.:-.     `.                             .:/                //
//                      \. /    `-._   `.""-----.,-..::(--"".\""`.  `:\                //
//                       `/         `-.__\          `-:\          `. `:\               //
//                                                                  `-._)              //
//                                                                                     //
//=====================================================================================//
```````

The SHARK library implements the soft contributions at low transverse momentum needed
for the NNLO expansion of the transverse momentum resummation formula, for processes
of the type Q+Qbar+F, where Q is a heavy quark and F a colour singlet system.

Full details and definitions can be found in the associated publication:

                 S. Devoto and J. Mazzitelli,
                 "Soft contributions to heavy quark production in arbitrary kinematics,"
                 [arXiv:2509.17509 [hep-ph]].
                 https://arxiv.org/abs/2509.17509
                     
Dependencies: 
       - GNU Scientific Library (GSL, https://www.gnu.org/software/gsl/gsl.html)
       - Polylogarithm (https://github.com/Expander/polylogarithm) [included in distribution]

Make sure the GSL library is installed and visible to your compiler (e.g. via `LD_LIBRARY_PATH` or system paths).

The code can be compiled from the main directory with the command "make".
The compiled library will be placed in the lib/ folder.
An example program is provided in the example/ directory.

Contact: javimazzitelli@gmail.com

