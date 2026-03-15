# NX 2306
# Journal created by andyc on Sun Mar 15 01:21:59 2026 Eastern Daylight Time
#
import math
import NXOpen
def main() : 

    theSession  = NXOpen.Session.GetSession()
    workPart = theSession.Parts.Work
    displayPart = theSession.Parts.Display
    # ----------------------------------------------
    #   Menu: Tools->Utilities->Expressions...
    # ----------------------------------------------
    markId1 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Visible, "Start")
    
    theSession.SetUndoMarkName(markId1, "Expressions Dialog")
    
    markId2 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Edit Expression")
    
    expression1 = workPart.Expressions.FindObject("casing_id")
    unit1 = workPart.UnitCollection.FindObject("MilliMeter")
    workPart.Expressions.EditExpressionWithUnits(expression1, unit1, "149.5")
    
    markId3 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Edit Expression")
    
    expression2 = workPart.Expressions.FindObject("casing_od")
    workPart.Expressions.EditExpressionWithUnits(expression2, unit1, "152.5")
    
    markId4 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Edit Expression")
    
    expression3 = workPart.Expressions.FindObject("chamber_length")
    workPart.Expressions.EditExpressionWithUnits(expression3, unit1, "149.8")
    
    markId5 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Edit Expression")
    
    expression4 = workPart.Expressions.FindObject("combustion_annulus_radial_height")
    workPart.Expressions.EditExpressionWithUnits(expression4, unit1, "20.8")
    
    markId6 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Edit Expression")
    
    expression5 = workPart.Expressions.FindObject("dil_length")
    workPart.Expressions.EditExpressionWithUnits(expression5, unit1, "59.9")
    
    markId7 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Edit Expression")
    
    expression6 = workPart.Expressions.FindObject("dilHole_qty")
    workPart.Expressions.EditExpressionWithUnits(expression6, NXOpen.Unit.Null, "33")
    
    markId8 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Edit Expression")
    
    expression7 = workPart.Expressions.FindObject("film_hole_area")
    unit2 = workPart.UnitCollection.FindObject("SquareMilliMeter")
    workPart.Expressions.EditExpressionWithUnits(expression7, unit2, "692")
    
    markId9 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Expressions")
    
    markId10 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Make Up to Date")
    
    objects1 = [NXOpen.NXObject.Null] * 7 
    objects1[0] = expression4
    objects1[1] = expression3
    objects1[2] = expression2
    objects1[3] = expression1
    objects1[4] = expression6
    objects1[5] = expression5
    objects1[6] = expression7
    theSession.UpdateManager.MakeUpToDate(objects1, markId10)
    
    markId11 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "NX update")
    
    nErrs1 = theSession.UpdateManager.DoUpdate(markId11)
    
    theSession.DeleteUndoMark(markId11, "NX update")
    
    theSession.DeleteUndoMark(markId10, None)
    
    theSession.DeleteUndoMark(markId9, None)
    
    theSession.SetUndoMarkName(markId1, "Expressions")
    
    markId12 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Visible, "Start")
    
    theSession.SetUndoMarkName(markId12, "Expressions Dialog")
    
    # ----------------------------------------------
    #   Dialog Begin Expressions
    # ----------------------------------------------
    markId13 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Edit Expression")
    
    expression8 = workPart.Expressions.FindObject("inner_dilHole_dia")
    workPart.Expressions.EditExpressionWithUnits(expression8, unit1, "6.95")
    
    markId14 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Edit Expression")
    
    expression9 = workPart.Expressions.FindObject("inner_film_hole_qty")
    workPart.Expressions.EditExpressionWithUnits(expression9, NXOpen.Unit.Null, "63")
    
    markId15 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Edit Expression")
    
    expression10 = workPart.Expressions.FindObject("inner_liner_annulus_gap")
    workPart.Expressions.EditExpressionWithUnits(expression10, unit1, "18.17")
    
    markId16 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Edit Expression")
    
    expression11 = workPart.Expressions.FindObject("inner_liner_id")
    workPart.Expressions.EditExpressionWithUnits(expression11, unit1, "78.3")
    
    markId17 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Edit Expression")
    
    expression12 = workPart.Expressions.FindObject("inner_liner_od")
    workPart.Expressions.EditExpressionWithUnits(expression12, unit1, "80.8")
    
    markId18 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Edit Expression")
    
    expression13 = workPart.Expressions.FindObject("inner_primHole_dia")
    workPart.Expressions.EditExpressionWithUnits(expression13, unit1, "4.67")
    
    markId19 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Edit Expression")
    
    expression14 = workPart.Expressions.FindObject("inner_secHole_dia")
    workPart.Expressions.EditExpressionWithUnits(expression14, unit1, "6.5")
    
    markId20 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Edit Expression")
    
    expression15 = workPart.Expressions.FindObject("outer_annulus_gap")
    workPart.Expressions.EditExpressionWithUnits(expression15, unit1, "11.9")
    
    markId21 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Edit Expression")
    
    expression16 = workPart.Expressions.FindObject("outer_dilHole_dia")
    workPart.Expressions.EditExpressionWithUnits(expression16, unit1, "8.54")
    
    markId22 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Edit Expression")
    
    expression17 = workPart.Expressions.FindObject("outer_film_hole_qty")
    workPart.Expressions.EditExpressionWithUnits(expression17, NXOpen.Unit.Null, "94")
    
    markId23 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Edit Expression")
    
    expression18 = workPart.Expressions.FindObject("outer_liner_id")
    workPart.Expressions.EditExpressionWithUnits(expression18, unit1, "122.62")
    
    markId24 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Edit Expression")
    
    expression19 = workPart.Expressions.FindObject("outer_liner_od")
    workPart.Expressions.EditExpressionWithUnits(expression19, unit1, "124.97")
    
    markId25 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Edit Expression")
    
    expression20 = workPart.Expressions.FindObject("outer_primHole_dia")
    workPart.Expressions.EditExpressionWithUnits(expression20, unit1, "5.23")
    
    markId26 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Edit Expression")
    
    expression21 = workPart.Expressions.FindObject("outer_secHole_dia")
    workPart.Expressions.EditExpressionWithUnits(expression21, unit1, "7.87")
    
    markId27 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Edit Expression")
    
    expression22 = workPart.Expressions.FindObject("prim_sec_length")
    workPart.Expressions.EditExpressionWithUnits(expression22, unit1, "45.0")
    
    markId28 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Edit Expression")
    
    expression23 = workPart.Expressions.FindObject("primHole_qty")
    workPart.Expressions.EditExpressionWithUnits(expression23, NXOpen.Unit.Null, "15")
    
    markId29 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Edit Expression")
    
    expression24 = workPart.Expressions.FindObject("secHole_qty")
    workPart.Expressions.EditExpressionWithUnits(expression24, NXOpen.Unit.Null, "15")
    
    markId30 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Edit Expression")
    
    expression25 = workPart.Expressions.FindObject("shaft_tunnel_id")
    workPart.Expressions.EditExpressionWithUnits(expression25, unit1, "39.0")
    
    markId31 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Edit Expression")
    
    expression26 = workPart.Expressions.FindObject("shaft_tunnel_od")
    workPart.Expressions.EditExpressionWithUnits(expression26, unit1, "42.0")
    
    markId32 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Edit Expression")
    
    expression27 = workPart.Expressions.FindObject("wall_thickness")
    workPart.Expressions.EditExpressionWithUnits(expression27, unit1, "1.6")
    
    markId33 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Expressions")
    
    markId34 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Make Up to Date")
    
    objects2 = [NXOpen.NXObject.Null] * 20 
    objects2[0] = expression26
    objects2[1] = expression25
    objects2[2] = expression27
    objects2[3] = expression19
    objects2[4] = expression18
    objects2[5] = expression15
    objects2[6] = expression12
    objects2[7] = expression11
    objects2[8] = expression10
    objects2[9] = expression23
    objects2[10] = expression20
    objects2[11] = expression24
    objects2[12] = expression13
    objects2[13] = expression14
    objects2[14] = expression21
    objects2[15] = expression8
    objects2[16] = expression16
    objects2[17] = expression22
    objects2[18] = expression9
    objects2[19] = expression17
    theSession.UpdateManager.MakeUpToDate(objects2, markId34)
    
    markId35 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "NX update")
    
    nErrs2 = theSession.UpdateManager.DoUpdate(markId35)
    
    theSession.DeleteUndoMark(markId35, "NX update")
    
    theSession.DeleteUndoMark(markId34, None)
    
    theSession.DeleteUndoMark(markId33, None)
    
    theSession.SetUndoMarkName(markId12, "Expressions")
    
    markId36 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Visible, "Start")
    
    theSession.SetUndoMarkName(markId36, "Expressions Dialog")
    
    # ----------------------------------------------
    #   Dialog Begin Expressions
    # ----------------------------------------------
    markId37 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Expressions")
    
    theSession.DeleteUndoMark(markId37, None)
    
    markId38 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Expressions")
    
    markId39 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "Make Up to Date")
    
    markId40 = theSession.SetUndoMark(NXOpen.Session.MarkVisibility.Invisible, "NX update")
    
    nErrs3 = theSession.UpdateManager.DoUpdate(markId40)
    
    theSession.DeleteUndoMark(markId40, "NX update")
    
    theSession.DeleteUndoMark(markId39, None)
    
    theSession.DeleteUndoMark(markId38, None)
    
    theSession.SetUndoMarkName(markId36, "Expressions")
    
    # ----------------------------------------------
    #   Menu: File->Save
    # ----------------------------------------------
    partSaveStatus1 = workPart.Save(NXOpen.BasePart.SaveComponents.TrueValue, NXOpen.BasePart.CloseAfterSave.FalseValue)
    
    partSaveStatus1.Dispose()
    # ----------------------------------------------
    #   Menu: Tools->Automation->Journal->Stop Recording
    # ----------------------------------------------
    
if __name__ == '__main__':
    main()