function out = isMRI(mri)

fieldList = [{'srcbext'   }
    {'analyzehdr'}
    {'bhdr'      }
    {'vol'       }
    {'niftihdr'  }
    {'fspec'     }
    {'pwd'       }
    {'flip_angle'}
    {'tr'        }
    {'te'        }
    {'ti'        }
    {'vox2ras0'  }
    {'volsize'   }
    {'height'    }
    {'width'     }
    {'depth'     }
    {'nframes'   }
    {'vox2ras'   }
    {'nvoxels'   }
    {'xsize'     }
    {'ysize'     }
    {'zsize'     }
    {'x_r'       }
    {'x_a'       }
    {'x_s'       }
    {'y_r'       }
    {'y_a'       }
    {'y_s'       }
    {'z_r'       }
    {'z_a'       }
    {'z_s'       }
    {'c_r'       }
    {'c_a'       }
    {'c_s'       }
    {'vox2ras1'  }
    {'Mdc'       }
    {'volres'    }
    {'tkrvox2ras'}];
out = all(isfield(mri,fieldList));