program godas_driver
  use godas_3dvar


  ! initialize
  call g3dv_init()

  call g3dv_run()
  
  call g3dv_final()

end program godas_driver
