      subroutine get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx)          772
      data sml0,eps0,big0,mnlam0,rsqmax0,pmin0,exmx0  /1.0e-5,1.0e-6,9.9    774
     *e35,5,0.999,1.0e-9,250.0/
      sml=sml0                                                              774
      eps=eps0                                                              774
      big=big0                                                              774
      mnlam=mnlam0                                                          774
      rsqmax=rsqmax0                                                        775
      pmin=pmin0                                                            775
      exmx=exmx0                                                            776
      return                                                                777
      entry chg_fract_dev(arg)                                              777
      sml0=arg                                                              777
      return                                                                778
      entry chg_dev_max(arg)                                                778
      rsqmax0=arg                                                           778
      return                                                                779
      entry chg_min_flmin(arg)                                              779
      eps0=arg                                                              779
      return                                                                780
      entry chg_big(arg)                                                    780
      big0=arg                                                              780
      return                                                                781
      entry chg_min_lambdas(irg)                                            781
      mnlam0=irg                                                            781
      return                                                                782
      entry chg_min_null_prob(arg)                                          782
      pmin0=arg                                                             782
      return                                                                783
      entry chg_max_exp(arg)                                                783
      exmx0=arg                                                             783
      return                                                                784
      end                                                                   785
