import textwrap

from peptools.chem import get_fasta_from_mol
from peptools.io.fasta import _is_input_fasta
from peptools.io.model import InputAttribute
from peptools.io.model import InputFactory
from rdkit import Chem


def test_is_fasta():
    """Verifies a string is interpreted as a FASTA string."""
    fasta = textwrap.dedent(
        """\
        >sp|P43220|GLP1R_HUMAN Glucagon-like peptide 1 receptor OS=Homo sapiens OX=9606 GN=GLP1R PE=1 SV=2
        MAGAPGPLRLALLLLGMVGRAGPRPQGATVSLWETVQKWREYRRQCQRSLTEDPPPATDL
        FCNRTFDEYACWPDGEPGSFVNVSCPWYLPWASSVPQGHVYRFCTAEGLWLQKDNSSLPW
        RDLSECEESKRGERSSPEEQLLFLYIIYTVGYALSFSALVIASAILLGFRHLHCTRNYIH
        LNLFASFILRALSVFIKDAALKWMYSTAAQQHQWDGLLSYQDSLSCRLVFLLMQYCVAAN
        YYWLLVEGVYLYTLLAFSVLSEQWIFRLYVSIGWGVPLLFVVPWGIVKYLYEDEGCWTRN
        SNMNYWLIIRLPILFAIGVNFLIFVRVICIVVSKLKANLMCKTDIKCRLAKSTLTLIPLL
        GTHEVIFAFVMDEHARGTLRFIKLFTELSFTSFQGLMVAILYCFVNNEVQLEFRKSWERW
        RLEHLHIQRDSSMKPLKCPTSSLSSGATAGSSMYTATCQASCS
    """
    )
    assert _is_input_fasta(fasta)


def test_not_fasta():
    """Verifies a string is not interpreted as a FASTA string."""
    assert not _is_input_fasta("CNCNCNCNC")


def test_input_factory():
    """Generates input data from a molecule."""
    input_data = InputFactory.new(Chem.MolFromSmiles("CC"), "my_mol")
    assert input_data[InputAttribute.MOL_NAME] == "my_mol"
    assert input_data[InputAttribute.MOL_OBJECT].GetNumAtoms() == 2
    assert input_data[InputAttribute.FASTA] == "X"
