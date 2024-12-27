test_that("proforma clean sequences works", {
    expect_identical(
        .proforma_clean_sequences(
            c("EM[+15.9949]EVEES[+79.9663]PEK", "EM[+15.995]EVEES[-18.01]PEK")
        ),
        c("EMEVEESPEK", "EMEVEESPEK")
    )
})

test_that("proforma delta masses without prefixes are supported", {
    expect_equal(
        .proforma_delta_masses(
            c("EMEVEESPEK", "EM[+15.995]EVEESPEK")
        ),
        list(
            rep(0, 10), 
            c(0, 15.995, 0, 0, 0, 0, 0, 0, 0, 0)
        )
    )
    expect_equal(
        .proforma_delta_masses(
            c("EM[+15.9949]EVEES[+79.9663]PEK", "EM[+15.995]EVEES[-18.01]PEK")
        ),
        list(
            c(0, 15.9949, 0, 0, 0, 0, 79.9663, 0, 0, 0),
            c(0, 15.995, 0, 0, 0, 0, -18.01, 0, 0, 0)
        )
    )
})

test_that("proforma delta masses with prefixes are supported", {
    expect_equal(
        .proforma_delta_masses(
            c("EM[U:+15.9949]EVEES[M:+79.9663]PEK",
              "EM[X:+15.995]EVEES[R:-18.01]PEK")
        ),
        list(
            c(0, 15.9949, 0, 0, 0, 0, 79.9663, 0, 0, 0),
            c(0, 15.995, 0, 0, 0, 0, -18.01, 0, 0, 0)
        )
    )
})
