def modify_fasta_according_to_cys_keys_pichemist(
    fasta_in,
    no_free_cys_thiols,
    n_disulfide_bonds_str,
):
    total_cys = fasta_in.count("C")

    if no_free_cys_thiols:
        n_cys_to_replace = total_cys

    else:
        if n_disulfide_bonds_str == "max":
            # Max disulfides means pairing as many cysteines as possible
            n_cys_to_replace = total_cys if total_cys % 2 == 0 else total_cys - 1

        else:
            try:
                n_disulfide_bonds = int(float(n_disulfide_bonds_str))
            except Exception:
                raise ValueError(
                    f"Cannot convert n_disulfide_bonds='{n_disulfide_bonds_str}' to integer."
                )

            if n_disulfide_bonds <= 0:
                n_cys_to_replace = 0
            else:
                n_cys_to_replace = n_disulfide_bonds * 2

                if n_cys_to_replace > total_cys:
                    raise ValueError(
                        f"Specified {n_disulfide_bonds} disulfide bonds "
                        f"but only {total_cys} cysteine residues are present."
                    )

    return fasta_in.replace("C", "X", n_cys_to_replace)


def modify_fasta_according_to_cys_keys_mec(
    fasta_in,
    no_free_cys_thiols,
    n_disulfide_bonds_str,
):
    total_cys = fasta_in.count("C")

    if n_disulfide_bonds_str == "max":
        # Pair as many cysteines as possible
        n_cys_in_disulfides = total_cys if total_cys % 2 == 0 else total_cys - 1
        n_disulfide_bonds = n_cys_in_disulfides // 2

    else:
        try:
            n_disulfide_bonds = int(float(n_disulfide_bonds_str))
        except Exception:
            raise ValueError(
                f"Cannot convert n_disulfide_bonds='{n_disulfide_bonds_str}' to integer."
            )

        if n_disulfide_bonds <= 0:
            n_cys_in_disulfides = 0
        else:
            n_cys_in_disulfides = n_disulfide_bonds * 2

    if n_cys_in_disulfides > total_cys:
        raise ValueError(
            f"Specified {n_disulfide_bonds} disulfide bonds "
            f"but only {total_cys} cysteine residues are present."
        )

    # Replace cysteines involved in disulfides
    fasta_out = fasta_in.replace("C", "𝒞", n_cys_in_disulfides)

    if no_free_cys_thiols:
        # Replace remaining free cysteines
        fasta_out = fasta_out.replace("C", "X")

    return fasta_out
