test_that("Get same dfs out applying lr_to_ra on output of ra_to_lr", {
   
   fixed_df <- data.frame("value" = 1/5,
                          "param" = "P_tilde",
                          k = 1,
                          j = 1:5)
   
   varying_df <- data.frame("value" = (1:5)/sum(1:5),
                            "param" = "P",
                            k = 1,
                            j = 1:5)
   
   expect_equal(lr_to_ra(fixed_df = fixed_df,
                         varying_lr_df =ra_to_lr(varying_df),
                         varying_df = varying_df),
                varying_df)
   
})

# test_that("Get same dfs out applying ra_to_lr on output of lr_to_ra",{
#    temp_varying <-
#       with(list_of_dfs,
#            lr_to_ra(fixed_df,varying_lr_df,varying_df))
#    temp_lr <- ra_to_lr(temp_varying)
#    
#    expect_equal(temp_lr,list_of_dfs$varying_lr_df)
# })
