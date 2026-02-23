# NX 2306
# Auto-generated Journal for Combustor Update
import math
import NXOpen

LinerID = "129.426"
LinerOD = "132.426"
Length = "218.503"
NumPrim = "6.0"
NumSec = "6.0"
NumDil = "6.0"
PrimDia = "8.279"
SecDia = "15.114"
DilDia = "19.012"
AnnulusGap = "21.187"

def main():
    theSession  = NXOpen.Session.GetSession()
    workPart = theSession.Parts.Work
    
    markId1 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Visible, "Start")
    unit1 = workPart.UnitCollection.FindObject("MilliMeter")
    
    # Update Liner ID
    expression1 = workPart.Expressions.FindObject("LID")
    workPart.Expressions.EditExpressionWithUnits(expression1, unit1, LinerID)
    
    # Update Liner OD
    expression2 = workPart.Expressions.FindObject("LOD")
    workPart.Expressions.EditExpressionWithUnits(expression2, unit1, LinerOD)
    
    # Update Hole Counts (Using NXOpen.Unit.Null for unitless)
    expression3 = workPart.Expressions.FindObject("numD")
    workPart.Expressions.EditExpressionWithUnits(expression3, NXOpen.Unit.Null, NumDil)
    
    expression4 = workPart.Expressions.FindObject("numP")
    workPart.Expressions.EditExpressionWithUnits(expression4, NXOpen.Unit.Null, NumPrim)
    
    expression5 = workPart.Expressions.FindObject("numS")
    workPart.Expressions.EditExpressionWithUnits(expression5, NXOpen.Unit.Null, NumSec)
    
    # Update Length
    expression6 = workPart.Expressions.FindObject("Length")
    workPart.Expressions.EditExpressionWithUnits(expression6, unit1, Length)
    
    # Update Dilution Dia
    expression7 = workPart.Expressions.FindObject("DDia")
    workPart.Expressions.EditExpressionWithUnits(expression7, unit1, DilDia)
    
    # Update Annulus Gap
    expression8 = workPart.Expressions.FindObject("AnnulusGap")
    workPart.Expressions.EditExpressionWithUnits(expression8, unit1, AnnulusGap)
    
    # Update Model
    markId12 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Make Up to Date")
    objects1 = [expression5, expression3, expression6, expression8, expression4, expression7, expression2, expression1]
    theSession.UpdateManager.MakeUpToDate(objects1, markId12)
    
    markId13 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "NX update")
    theSession.UpdateManager.DoUpdate(markId13)
    theSession.DeleteUndoMark(markId13, "NX update")

if __name__ == '__main__':
    main()
