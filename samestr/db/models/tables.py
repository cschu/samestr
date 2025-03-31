# pylint: disable=R0903

""" module docstring """

from sqlalchemy import Column, Integer, String, ForeignKey, Float
from sqlalchemy.orm import relationship

from .meta import Base


class Marker(Base):
    __tablename__ = "marker"
    
    id = Column(Integer, primary_key=True,)
    contig_id = Column(String, index=True,)
    gene_id = Column(String,)
    clade_id = Column(Integer, ForeignKey("clade.id"))
    begin = Column(Integer,)
    end = Column(Integer,)
    strand = Column(String,)
    seq = Column(String,)
    length = Column(Integer,)
    total_marker_length = Column(Integer,)
    
    clade = relationship("Clade", back_populates="markers")

    # contig_id, gene_id, beg, end, strand, seq.seq

class Clade(Base):
    __tablename__ = "clade"

    id = Column(Integer, primary_key=True,)
    name = Column(String, index=True,)

    markers = relationship("Marker", back_populates="clade")
